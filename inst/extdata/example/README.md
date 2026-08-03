# Minimal example data

`ModuleSelection/` contains a synthetic 10-node module for exercising
`build_module_manifest()`. It includes one self-loop and one `M0` edge so the
adapter's exclusion and QC behavior can be verified without patient-level or
study-derived data.

After installing the R package:

```r
example_root <- system.file("extdata", "example", "ModuleSelection", package = "MLSnpDR")
output_dir <- tempfile("mlsnpdr-manifest-")
manifest <- build_module_manifest(example_root, output_dir)
manifest
```

Future releases will extend this directory with reduced synthetic inputs for
the remaining pipeline stages.
