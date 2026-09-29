# FI×sleep LAVA resource and throughput telemetry — 2026-09-27 01:28 UTC

Read-only host/process samples; no worker controls or analysis inputs were changed.

## Active execution

At 01:28:17 UTC the coordinator reported 18,833 receipts; at 01:28:29 it reported 18,835. Six workers are requested and two are effective under the existing adaptive fallback. At the second sample, worker 72534's live R child was executing `insomnia` locus 1809 and worker 72535's was executing `snoring` locus 1808; both matched separate atomic claims. R child PIDs were 64527 and 64520, respectively. No launch failures or stale recoveries were recorded.

The 12-second window added two receipts (nominally 600/hour if extrapolated); this is too short and noisy to use for an ETA.

## Host resource samples

| UTC | Swap used / total | System memory free | CPU | Physical memory | R children |
|---|---:|---:|---|---|---|
| 01:28:17 | 3,574 / 5,120 MiB | 60% | 97.66% user, 2.10% system, 0.23% idle | 7,534 MiB used, 97 MiB unused | 269.4% and 265.3% CPU; 692,176 and 759,008 KiB RSS |
| 01:28:29 | 3,470 / 5,120 MiB | 60% | 63.73% user, 15.41% system, 20.85% idle | 7,512 MiB used, 119 MiB unused | 100.0% and 98.8% CPU; 532,768 and 483,696 KiB RSS |

Load averages at the second sample: 11.00, 12.63, 21.02 on the shared host. A brief disk sample at 01:23 UTC showed external volume `disk4` at 37.61 MiB/s; the external SSD had 1.5 TiB free and the internal volume 3.8 GiB free. Disk throughput is a point sample, not an average.

## Measured throughput context

A previous serial-era short-window audit measured approximately 81 receipts/hour. The 00:15–00:45 full audit interval at two effective parallel workers measured approximately 310 validated receipts/hour, about 3.8 times that serial-era rate. Both are observed short-window rates on a shared host, not controlled benchmarks. At the latest full audit, 11,205 slots remained; a roughly 310/hour continuation would imply about 36 hours, strictly provisional. The current 12-second count change is not substituted for that longer measurement.

## Safety and scientific status

The adaptive coordinator remains at two effective workers. Although swap briefly dipped below its absolute downshift threshold in the second sample, the machine had only 119 MiB physically unused, high CPU activity, and substantial existing swap use. Do not raise concurrency based on this short fluctuation. Six total workers have not been sustained. The lock and prepared-input hashes remain unchanged. All 12 locked LAVA 1% trait gates fail, so no local-sharing inference is admissible.
