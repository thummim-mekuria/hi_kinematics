# Function variants

One file per duplicated function. Each contains every definition found in the
notebook exports, in notebook order, with a header saying where it came from
and whether it is a new version or a repeat.

| function | defs | notebooks | distinct bodies | distinct exact texts |
|---|---:|---:|---:|---:|
| [`fit_gradient`](fit_gradient.py) | 18 | 16 | 11 | 12 |
| [`moment_maps`](moment_maps.py) | 8 | 8 | 6 | 6 |
| [`make_zxz`](make_zxz.py) | 7 | 7 | 5 | 5 |
| [`read_DHT`](read_DHT.py) | 12 | 12 | 4 | 4 |
| [`plot_rvsv`](plot_rvsv.py) | 5 | 5 | 4 | 4 |
| [`first_moment`](first_moment.py) | 9 | 9 | 3 | 4 |
| [`surface_density`](surface_density.py) | 5 | 5 | 3 | 3 |
| [`plane_rvld`](plane_rvld.py) | 17 | 16 | 2 | 2 |
| [`get_subcube`](get_subcube.py) | 9 | 9 | 1 | 1 |
| [`get_theta`](get_theta.py) | 8 | 8 | 1 | 1 |
| [`pearson`](pearson.py) | 8 | 8 | 1 | 1 |
| [`second_moment`](second_moment.py) | 8 | 8 | 1 | 1 |
| [`separate_clouds`](separate_clouds.py) | 7 | 7 | 1 | 1 |

`distinct bodies` ignores indentation, blank lines and comment-only lines.
`distinct exact texts` is whitespace-sensitive. Where the two differ, copies
vary only in formatting or comments.

- **first_moment** — 148, 152, 154, 156, 157, 158, 159, 164, 165
- **fit_gradient** — 132, 134, 135, 137, 138, 141, 148, 152, 154, 156, 157, 158, 159, 162, 164, 165
- **get_subcube** — 148, 152, 154, 156, 157, 158, 159, 164, 165
- **get_theta** — 148, 154, 156, 157, 158, 159, 164, 165
- **make_zxz** — 152, 156, 157, 158, 159, 164, 165
- **moment_maps** — 148, 154, 156, 157, 158, 159, 164, 165
- **pearson** — 148, 152, 154, 155, 156, 158, 162, 165
- **plane_rvld** — 132, 134, 135, 137, 138, 141, 148, 152, 154, 156, 157, 158, 159, 162, 164, 165
- **plot_rvsv** — 152, 154, 158, 162, 165
- **read_DHT** — 027, 028, 029, 030, 031, 033, 034, 035, 036, 037, 038, 039
- **second_moment** — 148, 154, 156, 157, 158, 159, 164, 165
- **separate_clouds** — 044, 045, 046, 047, 048, 049, 050
- **surface_density** — 022, 023, 024, 041, 045
