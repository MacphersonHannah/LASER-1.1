# calibrated_thetas_reference.csv

Reference table of previously calibrated theta values, for `theta_mode = 'lift'`.
You can use the provided thetas, or add a row for a new tissue/dataset and use it immediately with no code changes needed.

## Columns

| column             | meaning                                                              |
|--------------------|-----------------------------------------------------------------------|
| `source`           | unique name you reference via `theta_lift_source` / `theta_lift_table` |
| `theta`            | the calibrated value                                                 |
| `achieved_fpr`     | empirical FPR the calibration converged on (should be near your nominal target, e.g. 0.05) |
| `n_loci`           | number of loci the calibration was run over                          |
| `n_samples_pooled` | number of samples pooled into that calibration run                   |
| `calibration_date` | when it was calibrated (YYYY-MM-DD)                                   |
| `laser_version`    | LASER version/tag used                                               |
| `notes`            | anything else worth recording (cohort, tissue detail, caveats)       |

## Adding your own

1. Run the pipeline once with `theta_mode = 'calibrate_pooled'` on the
   sample(s) you want to calibrate on.
2. Read the resulting theta and achieved FPR from
   `allelic_effect_testing/_theta/theta_pooled.json`.
3. Append a new row to this CSV with a `source` name that isn't already
   used (e.g. `liver_v1`, `my_lab_kidney`).
4. Reference it: `--theta_lift_source liver_v1`, or per-sample via a
   `theta_lift_table` CSV (`sample,source`).

`source` names must be unique - the pipeline will error out with a clear
message if it finds a duplicate, rather than silently picking one.
