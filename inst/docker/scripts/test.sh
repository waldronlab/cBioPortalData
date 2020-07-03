#/usr/bin/env bash
Rscript -e 'res=devtools::test(reporter="summary");df=as.data.frame(res);if(sum(df$failed) > 0 || any(df$error)) {q(status=1)}'
