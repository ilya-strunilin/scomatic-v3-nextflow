# Compute1 Nextflow + LSF configuration

`nextflow.config` defines the scientific defaults and the `compute1` profile
adds the LSF executor. The profile's public default is:

```groovy
process.executor = 'lsf'
process.queue = 'general'
workDir = "<run_root>/work"
```

Every task writes its LSF combined log to
`<run_root>/logs/nextflow_task_<job-id>.log`; the controller has a separate
`nextflow_controller_<job-id>.log`. A real `run_root` is required specifically
to avoid `null/logs/...` paths and associated mail notifications.

The task container is configurable through `lsf_container`. The default is the
image used by this workflow's Compute1 deployment; change it if your account
uses a different approved image. Similarly, set `lsf_group` only if your LSF
account requires a group, and export `LSF_DOCKER_VOLUMES` before submission if
the Compute1 Docker integration requires explicit bind mounts for your data,
reference, scratch, and repository paths.

The controller script sets `NXF_HOME` and `NXF_WORK` inside `run_root`, so
Nextflow state, launch cache, and controller work do not accumulate in HOME.
It does not bundle or publish any LSF executables, certificates, local queue
files, credentials, or site configuration.

Use `--resume` only with the identical run root and inputs. The task wrappers
also recognize verified output markers and reuse completed donor stages rather
than re-running them.
