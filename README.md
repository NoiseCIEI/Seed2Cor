# Seed2Cor

Parallel single component ambient noise correlation from `SEED`
- Extract `SAC` & response from `SEED`
- Remove instrument response
- Temporal (1-bit, run-absolute-mean, cut events) &
  spectral (whiten) normalization
- Cross/Auto-correlate & stack

# Directory structures

- `input/`: template files for input
- `src/`: source code
- `tools/`: auxiliary scripts

# Usage

## Dependencies

Currently tested on

- RHEL (Red Hat Enterprise Linux): 6.9, 7.3
- gcc: 6.1.0, 4.4.7 20120313 (Red Hat 4.4.7-18)
- FFTW: 3.3.7, 3.3.4 (failed!), 3.2.3
- rdseed 5.3.1
- evalresp 3.3.3

## Compile

```
$ cd src/
$ make
```

The output could be compared with `src/make.log`.

## Input

- `parameters`: processing parameters
  - see `input/parameters.lst` for template
- `seed.lst`: list of `SEED` files
  - see `input/seed.lst` for template
  - use `tools/generate_seedlst.sh` to generate
- `station.lst`: list of stations
  - see `input/station.lst` for template
  - use `tools/generate_stationlst.sh` to generate
  - optional flag
    ```
    0-N
    For stations with same number n (except 0), do cc,
    also do cc with group 0,
    but do not do cc with other numbers.
    For number 0, do cc with all other group,
    but do not do cc with all other group 0 members.
    ```

## Run

```shell
$ Seed2Cor parameter_file [nthreads]
```

A `Python` wrapper `tools/s2c_wrapper.py` with a `YAML` parameter file
`input/param.yml` might be helpful.

## Output

### File Structure

- `yyyy.MMM/`
  - `yyyy.MMM.dd/`
  - `COR/`: monthly correlations
    - `STA/`
      - `COR_MASTER_SLAVE.SAC`
  - `COR_D/`: daily correlations
    - `STA/`
      - `COR_MASTER_SLAVE_dd.SAC`
  - `Cor_dayflag.lst`: days if data available

### SAC Header

- `KEVNM`: source station
- `KSTNM`: receiver station
- `USER0`: # of days stacked

## Attention

- The temperal normalization method of 'earthquake cutting' has NOT been well tested. Use with caution!
- The spectrum reshape function is currently removed from the code. (`Smooth2()` in `Whiten.c` does nothing)

# Reference

Bensen, G. D., Ritzwoller, M. H., Barmin, M. P., Levshin, A. L., Lin, F., Moschetti, M. P., et al. (2007). Processing seismic ambient noise data to obtain reliable broad-band surface wave dispersion measurements. *Geophysical Journal International*, 169(3), 1239–1260.

Lin, F.-C., Ritzwoller, M. H., & Shen, W. (2011). On the reliability of attenuation measurements from ambient noise cross-correlations. *Geophysical Research Letters*, 38(11).
