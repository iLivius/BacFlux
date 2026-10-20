# Contributing to BacFlux

Contributions are welcome. Three rules apply to every change.

## Don't bundle restricted databases

Never commit a third-party database to this repository. Put its download URL in
`config/config.yaml` instead, and say in the documentation which licence it
comes under. Each user then downloads it under the provider's own terms.

## Don't copy non-commercial code

BacFlux is MIT-licensed. Code released under CC BY-NC-SA, or any other
non-commercial or share-alike licence, must not be copied into it. That holds
even for a part that is switched off by default: the licence travels with the
code, whether or not the code runs. Reading such code to understand how
something works, and then writing your own version, is fine.

[Licensing](https://ilivius.github.io/BacFlux/about/licensing/) explains how
BacFlux keeps its own code apart from third-party code and data.

## Leave an audit trail for every filtering decision

Whenever a step removes, keeps or reclassifies something (a contig, a gene
call, a predicted mobile element), it writes a tab-separated file with one row
per decision and a column giving the reason. Nothing is dropped silently.

## Before you open a pull request

Run the tests described under
[For contributors](https://ilivius.github.io/BacFlux/reference/running/#for-contributors).
