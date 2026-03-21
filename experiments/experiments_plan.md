Сейчас важнее всего не “найти лучший конфиг”, а правильно поставить поиск. У тебя пока нет gold standard, значит в первой итерации надо оптимизировать не один скалярный score, а искать `Pareto front`: много находок, но без взрыва мусора. И я бы точно не делал 50 полных прогонов на все 500+ генов. Правильнее так: `tuning panel` из 60-100 генов, потом 3-5 лучших конфигов валидировать на всем наборе. Иначе ты переедешь по compute, а выводы будут хуже.

Важный момент по порядку тюнинга: `ortholog scope` и `read_len/step` надо выбирать раньше BAM thresholds и раньше окончательной настройки VarScan. Причина простая: они меняют сам состав ортологов и покрытие псевдоридов, а значит меняют распределение сигналов, на которых потом работают BAM-фильтры и caller. То есть если сначала тюнить `min_mapped_pct_of_generated` или `min_var_freq`, а потом менять `scope` или `step`, ты фактически заново меняешь задачу.

**Первый ран**
Я бы не начинал с `25/25`. Для BWA это уже слишком коротко для такой задачи: резко вырастет неоднозначное выравнивание, особенно на консервативных/повторных участках и при кросс-видовых ортологах. Перекрытие ридов полезно, и не только для “больше покрытия”: оно дает несколько сдвинутых контекстов вокруг одной позиции, что помогает и выравниванию, и pileup. Поэтому первый recall-run я бы делал так:

- `ortholog_selection.scope = all`
- `read_generation.read_len = 50`
- `read_generation.step = 25`
- `read_generation.pseudo_read_phred = 30`
- `bam_filtering.enabled = false`
- `variant_calling.min_coverage = 1`
- `variant_calling.min_reads2 = 1`
- `variant_calling.min_var_freq = 0.01`

То есть максимально либеральный caller, но без ухода в совсем бессмысленный `25bp`. И параллельно я бы почти сразу запускал второй anchor-run: тот же caller, но `75/35`. Тогда ты сразу увидишь, где основной источник “новых” мутаций: более короткие риды или просто открытый VarScan.

**Что считать по каждому run**
Одного score сейчас не надо. Я бы для каждого run собирал минимум:

1. `success_rate` по генам.
2. `genes_with_any_variant`.
3. `median variants/gene` и `p95 variants/gene`.
4. `fraction zero-variant genes`.
5. `ClinVar-overlap count` и отдельно `PLP/LP overlap`, если потом распарсишь `CLNSIG`.
6. `gnomAD-annotated fraction`.
7. `ortholog source mix` из `ortholog_resolution.csv` (`ncbi_scope`, `ncbi_all`, `blast`).
8. Если BAM filtering включен: `dropped_homologues`, `pct_kept`, `mapped_pct_of_generated_total`.

Это даст тебе не “лучший run”, а картину: recall, шум, устойчивость.

**План на 47 запусков**
Ниже не grid, а staged search. Он учитывает взаимодействия там, где они реально сильны, и не тратит слоты на слабые факторы вроде `pseudo_read_phred`.

1. Подготовка panel.
Возьми `80` генов для тюнинга и `100-150` в holdout. Tuning panel лучше стратифицировать по уже имеющемуся baseline:
`0 variants`, `low`, `medium`, `very high`, плюс разные источники ортологов из `ortholog_resolution.csv`.

2. Фаза A, `6` anchor-runs.
Нужны, чтобы быстро понять форму пространства.
`A1` текущий baseline.
`A2` ultra-recall: `all`, `50/25`, BAM off, `1/1/0.01`.
`A3` ultra-recall: `all`, `75/35`, BAM off, `1/1/0.01`.
`A4` ultra-recall: `mammals`, `50/25`, BAM off, `1/1/0.01`.
`A5` diagnostic: `eukaryota`, `50/25`, BAM off, `1/1/0.01`.
`A6` structural-only BAM: `all`, `75/35`, BAM on with `wrong_strand=true, lis=true, overlap=true`, thresholds permissive `0/100/0`, caller still `1/1/0.01`.

3. Фаза B, `10` run на `scope/read geometry`.
Здесь BAM off, caller фиксирован `1/1/0.01`.
`B1-B4`: `all` с `(50,25)`, `(75,25)`, `(75,35)`, `(100,50)`.
`B5-B8`: `mammals` с теми же четырьмя парами.
`B9-B10`: diagnostic no-overlap `(75,75)` для `all` и `mammals`.
После этой фазы выбираешь 1 лучший `scope` и 1-2 лучшие пары `read_len/step`.
Если `eukaryota` из `A5` слаб или почти копирует `all`, дальше его не трогаешь.

4. Фаза C, `7` run на структурные BAM-фильтры.
Берешь лучший `scope` и лучшую геометрию из фазы B. Thresholds пока держишь `0/100/0`, caller все еще либеральный.
Пробуешь все meaningful combinations stage toggles:
`wrong_strand only`
`lis only`
`overlap only`
`wrong_strand + lis`
`wrong_strand + overlap`
`lis + overlap`
`wrong_strand + lis + overlap`
Конфиг “BAM off” уже у тебя есть из прошлой фазы, его переиспользуешь как reference.

5. Фаза D, `7` run на VarScan.
Фиксируешь лучший `scope`, лучшую геометрию и лучший structural BAM mode. Thresholds BAM пока еще permissive.
Прогоняешь такую лестницу:
`(1,1,0.01)`
`(2,1,0.01)`
`(4,1,0.02)`
`(4,2,0.05)`
`(6,2,0.05)`
`(8,2,0.10)`
`(10,3,0.20)`
Здесь тройка = `(min_coverage, min_reads2, min_var_freq)`.
Это не полная сетка, но хорошо показывает, где начинается потеря полезного сигнала.

6. Фаза E, `9` run на BAM thresholds с interaction.
Фиксируешь уже выбранные `scope`, geometry, structural BAM combo и caller.
Запускаешь 9-run Latin-like design по:
`(min_mapped_pct_of_generated, max_pct_filtered, min_kept_pct_of_reference)`
с точками:
`(0,100,0)`
`(0,75,10)`
`(0,50,20)`
`(10,100,10)`
`(10,75,20)`
`(10,50,0)`
`(20,100,20)`
`(20,75,0)`
`(20,50,10)`
Это лучше, чем по одному крутить каждый порог: interaction между ними тут реально есть.

7. Фаза F, `3` local refinement.
Берешь 2-3 лучших конфига с Pareto front и делаешь вокруг них маленький ручной локальный search:
один run чуть мягче caller,
один чуть жестче BAM,
один компромиссный.
Это закрывает interaction `caller x BAM`, которую мы не добивали полным факториалом.

8. Фаза G, `5` validation-run на всех 500+ генах.
`G1` текущий baseline.
`G2` лучший recall-oriented конфиг.
`G3` лучший balanced конфиг.
`G4` лучший conservative конфиг.
`G5` либо второй balanced candidate, либо переключение `all <-> mammals` для winner-конфига.

Итого: `6 + 10 + 7 + 7 + 9 + 3 + 5 = 47` запусков.

**Чего бы я не делал в первой кампании**
1. Не тратил бы слоты на `pseudo_read_phred`. Пока держать `30`.
2. Не лез бы в `25/25`, если только один diagnostic run потом, когда станет ясно, что короткие риды реально помогают.
3. Не крутил бы одновременно `scope`, `geometry`, `BAM thresholds` и `VarScan`. Это даст не оптимизацию, а кашу.
4. Не выбирал бы “лучший конфиг” только по числу найденных вариантов. Нужен минимум 4-мерный взгляд: `recall proxy`, `noise proxy`, `run success`, `source stability`.

Если хочешь, следующим сообщением я могу превратить это в очень конкретную таблицу из `47` запусков в формате `run_id -> config diff from baseline`, чтобы тебе было удобно реально запускать.