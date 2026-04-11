## Test environments
* local Windows 11 install, R 4.4.1
* win-builder (devel)

## R CMD check results
0 errors | 0 warnings | 1 note

* NOTE: checking CRAN incoming feasibility
  * This is a new submission.

## URL checks
* NOTE: `urlchecker::url_check()` reports a 403 Forbidden error for 
`https://www.decathlon.fr/`. This is a false positive. The website uses 
strict anti-bot protections (like Cloudflare) that block automated HTTP HEAD 
requests from R. 
This URL is perfectly valid and accessible via a standard web browser.
