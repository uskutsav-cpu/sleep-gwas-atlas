# LAVA host resource telemetry — 2026-09-26 09:28 UTC

One read-only host sample was collected under the user's explicit LAVA monitoring instruction after sandbox process-inspection commands were denied. No process was stopped and no runner state was changed.

- **CPU:** 8 cores; load average 6.51/6.19/6.06; sampled CPU use 60.84% user, 23.11% system, 16.3% idle.
- **Memory:** 8 GiB total; `memory_pressure` reported 50% free; `top` showed 7,481 MiB used, 152 MiB unused, 1,435 MiB wired, 2,164 MiB compressed.
- **Swap:** 1,519.31/3,072 MiB used (1,552.69 MiB free), about 17 MiB below the scheduler's 1,536 MiB scale-down rung.
- **Active LAVA R children:** worker 2's `longsleep` job used 108% CPU and 393,360 KiB RSS; worker 1's `sleep_efficiency` job used 5% CPU and 130,432 KiB RSS. These were the assignments at sample time, not the later current claims.
- **Disk I/O:** the initial `iostat` row was 100.36 MB/s internal and 32.23 MB/s external; the next two one-second intervals were 0.26/0.40 and 0.08/0.27 MB/s.
- **Free disk:** 6.5 GiB internal, 1.6 TiB external.

The host had limited CPU idle, very little unused physical RAM, and swap just below the configured scale-down threshold. The current effective two-worker setting was preserved; no safe basis was established for increasing concurrency.

## Follow-up sample — 2026-09-26 09:41 UTC

- **CPU:** 8 cores; load averages 6.21/6.33/6.22; sampled CPU use 37.93% user, 13.79% system, 48.27% idle.
- **Memory:** `memory_pressure` reported 53% free; `top` showed 7,412 MiB used, 218 MiB unused, 1,397 MiB wired, and 2,146 MiB compressed.
- **Swap:** 1,519.31/3,072 MiB used (1,552.69 MiB free), unchanged from 09:28 and 16.7 MiB below the scheduler's 1,536 MiB scale-down rung.
- **Active LAVA R children:** worker 2's `sleepdur` locus 1096 used 549.2% CPU and 914,576 KiB RSS; worker 1's `sleep_apnea` locus 667 used 130.7% CPU and 436,928 KiB RSS.
- **Disk I/O:** the initial cumulative row was 99.54 MB/s internal and 31.95 MB/s external; the next one-second interval was 0.06 and 0.09 MB/s respectively.
- **Free disk:** 6.5 GiB internal, 1.6 TiB external.

At this sample, the two active R processes occupied about 1.3 GiB RSS. The existing two-worker configuration was retained because the recorded swap use was close to the scale-down threshold and the campaign's previous six-worker trial had already crossed the safeguards.

## Follow-up sample — 2026-09-26 09:49 UTC

- **CPU:** 8 cores; load averages 5.21/6.02/6.11; sampled CPU use 27.50% user, 15.62% system, 56.87% idle.
- **Memory:** `memory_pressure` reported 52% free; `top` showed 6,615 MiB used, 1,016 MiB unused, 1,427 MiB wired, and 2,142 MiB compressed.
- **Swap:** 1,511.31/3,072 MiB used (1,560.69 MiB free), about 25 MiB below the 1,536 MiB downshift rung.
- **Active LAVA R children:** worker 1's `snoring` locus 671 used 104.8% CPU and 223,488 KiB RSS; worker 2's `sleepiness` locus 1101 used 24.4% CPU and 130,832 KiB RSS.
- **Disk I/O:** initial cumulative rates were 99.03 MB/s internal and 31.78 MB/s external; the next one-second interval was 43.31 and 7.83 MB/s respectively.
- **Free disk:** most recent measured 6.5 GiB internal and 1.6 TiB external.

Memory pressure was lower than at 09:41, but swap remained near the scheduler's threshold. The effective two-worker setting was preserved; current measured throughput was already ahead of earlier four-worker intervals.
