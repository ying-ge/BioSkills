# CHANGES IN VERSION 2.18.0

## NEW FUNCTIONS
- `pathways` and `plotPathways`for summarizing & visualizing pathways Issue: [956](https://github.com/PoisonAlien/maftools/issues/956)
- `coGisticChromPlot` for plotting two GISTIC objects side-by-side. PR by [biosunsci](https://github.com/biosunsci) [954](https://github.com/PoisonAlien/maftools/pull/954)
- `readGistic` can take gistic output directory as an input. PR by [biosunsci](https://github.com/biosunsci) [954](https://github.com/PoisonAlien/maftools/pull/954)

## BUG FIXES
- Bug fixes while processing custom pathways
- Bug fix in `oncoplot` for drawing borders. [958](https://github.com/PoisonAlien/maftools/issues/958)
- Bug fix in `plotSignatures` for hardcoded axis limits. [949](https://github.com/PoisonAlien/maftools/issues/949)
- Bug fix in `mafSurvival` legend when samples argument is give. [937](https://github.com/PoisonAlien/maftools/issues/937)
- Bug fix in `subsetMaf` while handling only CNV events. [908](https://github.com/PoisonAlien/maftools/issues/908)
- Error handling when no deep/shallow CNV events found. [899](https://github.com/PoisonAlien/maftools/issues/899)
- Bug fix in `oncoplot` for duplicated values in gene list. [889](https://github.com/PoisonAlien/maftools/issues/889)

## ENHANCEMENTS
- Added argument `collapsePathway` to `oncoplot`. Issue: [956](https://github.com/PoisonAlien/maftools/issues/956)
- Improved `annovarToMaf` with better handling of indels and `Variant_Type`. Issue: [940](https://github.com/PoisonAlien/maftools/issues/940)
- Include absolute contribution of each signature in `extractSignatures` output. Issue: [939](https://github.com/PoisonAlien/maftools/issues/939)
- Added `tsbToPIDs` for custom names in `oncoplot`. Issue: Issue: [936](https://github.com/PoisonAlien/maftools/issues/936)
- Added `DSEL` protein to the database. Issue: [933](https://github.com/PoisonAlien/maftools/issues/933)
- Added `MUC3A` protein to the database. Issue: [932](https://github.com/PoisonAlien/maftools/issues/932)
- Added `showOnlyPathway` argument to `oncoplot`
- Added `pathdb` argument to `PlotOncogenicPathways`. Issue: [923](https://github.com/PoisonAlien/maftools/issues/923)
- Emit warnings when fishers test can not be performed during `somaticInteractions`. Issue: [921](https://github.com/PoisonAlien/maftools/issues/921)
- Added `leftMar` and `topMar` arguments to `somaticInteractions`. Issue: [913](https://github.com/PoisonAlien/maftools/issues/913)
- Added `toptBarLims` argument to oncoplot. Issue: [910](https://github.com/PoisonAlien/maftools/issues/910)
- Added `data` argument to `lollipopPlot` function. Issue: [894](https://github.com/PoisonAlien/maftools/issues/894)
- Added `sortByM1` and `sortByM2` argument to `coOncoplot`. Issue: [888](https://github.com/PoisonAlien/maftools/issues/888)
- Added arguments `leftBarVline`, `leftBarVlineCol`, `rightBarVline`, `rightBarVlineCol` `topBarHline` `topBarHlineCol` to `oncoplot`. Issue: [874](https://github.com/PoisonAlien/maftools/issues/874)
- Added `revPal` argument to `somaticInteractions`. Issue: [859](https://github.com/PoisonAlien/maftools/issues/859)
- Fix legend and color codes for numeric annotations in `oncoplot`. Issue: [363](https://github.com/PoisonAlien/maftools/issues/363)
