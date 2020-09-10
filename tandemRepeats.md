---
title: "Tandem repeats"
output: 
  html_document: 
    keep_md: yes
editor_options: 
  chunk_output_type: console
---


```r
requireNamespace("rtracklayer")
```

```
## Loading required namespace: rtracklayer
```

```r
suppressPackageStartupMessages(library("GenomicRanges"))
```



```r
loadAndReportChecksum <- function (file) {
  cat(system2("md5sum", file, stdout = TRUE))
  rtracklayer::import.gff(file)
}

tantan100  <- loadAndReportChecksum("~/FromSango/I69-4.tantan.f4.w100.gff3")
```

```
## 03b68db67012aa310061f8d98bb7664f  /home/charles/FromSango/I69-4.tantan.f4.w100.gff3
```

```r
tantan2000 <- loadAndReportChecksum("~/FromSango/I69-4.tantan.f4.w2000.gff3")
```

```
## 3fea3aa61955073d719357c0c58def60  /home/charles/FromSango/I69-4.tantan.f4.w2000.gff3
```

```r
ultra100   <- loadAndReportChecksum("~/FromSango/I69-4.ultra.p100.gff3")
```

```
## 864535e7e8bfb34a9290825b8dae2faa  /home/charles/FromSango/I69-4.ultra.p100.gff3
```

```r
ultra2000  <- loadAndReportChecksum("~/FromSango/I69-4.ultra.p2000.gff3")
```

```
## 846499539d0c3e2902b32f8c544a1a9b  /home/charles/FromSango/I69-4.ultra.p2000.gff3
```


```r
totLength <- function(X) sum(width(reduce(X)))
sapply(list(tantan100, tantan2000, ultra100, ultra2000), totLength)
```

```
## [1] 1459315 1634739 2334241 2650251
```

```r
totLength(intersect(tantan100, tantan2000))
```

```
## [1] 1383642
```

```r
totLength(intersect(ultra100, ultra2000))
```

```
## [1] 2327726
```

```r
totLength(setdiff(tantan2000, tantan100))
```

```
## [1] 251097
```
