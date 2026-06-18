# Changelog

## v1.2.0

### Changes

- Added plugin self-tests for cross-section reference values.
- Added plugin self-tests for sampled scattering events.
- Added plugin self-tests for spherical Bessel `j_n(x)` values.
- Added plugin self-tests for modified Bessel `K1` values.
- Enabled plugin self-tests in GitHub Actions CI.

## v1.1.0

### Changes

- Added spherical Bessel `j_n(x)` support for the Young-Koppel model.
- Added modified Bessel `K1` support for the Egelstaff diffusion model.
- Added GitHub Actions CI.

### Notes

- Improves compiler portability by removing reliance on C++ standard-library special Bessel functions.

## v1.0.0

- Initial release.
