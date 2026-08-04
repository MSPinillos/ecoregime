# Version 0.4.0

## Test environments

* Local:
  - macOS Tahoe 26.5, R 4.6.0 (x86_64-apple-darwin20)
* win-builder:
  - Windows Server 2022 x64, R Under development (2026-06-06 r90114 ucrt)(x86_64-w64-mingw32)
* mac-builder (https://mac.r-project.org/macbuilder/submit.html):
  - macOS Tahoe 26.6, R 4.6.1 Patched (aarch64-apple-darwin23)
* GitHub Actions:
  - macOS Tahoe 26.4, R 4.6.1 (aarch64-apple-darwin23)
  - Windows Server 2022 x64, R 4.6.1 (x86_64-w64-mingw32)
  - Ubuntu 24.04.4 LTS, R Under development (2026-06-21 r90185) (x86_64-pc-linux-gnu)
  - Ubuntu 24.04.4 LTS, R 4.6.1 (x86_64-pc-linux-gnu)
  - Ubuntu 24.04.4 LTS, R 4.5.3 (x86_64-pc-linux-gnu)

## R CMD check results

0 errors | 0 warnings | 1 notes

* checking examples ... [54s] NOTE
Examples with CPU (user + system) or elapsed time > 10s
                    user system elapsed
resilience_metrics 20.41   0.72   21.18
dist_edr            9.40   0.64   10.06

resilience_metrics includes examples for four functions.

## Reverse dependencies

* There are currently no downstream dependencies for this package.
