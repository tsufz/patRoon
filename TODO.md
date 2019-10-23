# Release

## general
- test negative subset indices
- refs to OpenBabel?
- convertMSFiles()
    - Agilent .d is also a directory?
    - Support OpenMS vendor conversion? (eg thermo)
- runWithoutCache? runWithCacheMode()? shortcut to withr::with_options(patRoon.cache.mode=...)

## AutoID

- ID level rules
    - GenForm scoring: somehow exclude non MS/MS candidates if MS/MS candidates are present?
        - already fine now with minimum rank?
        - make OR'ing possible?
    - add scorings for SIRIUS/DA
    - remove mustExist?
    - convert to YAML
        - no numeric level handling? makes file easier.
        - Filters could take perhaps order of file, eg if bestHit="4a" then ignore all rules below 4a
            - but this needs the rules...
            - maybe force a format? could be done while loading file, eg "[0-9]+[a-z]?"
- interface
    - groupFeaturesScreening() --> screenSuspects()?
    - deprecate old screen functions
    - also convert TASQ?
    - newProject(): update for new interface
    - annotateSuspects() --> annotate() latter is a function (but not generic) from ggplot2 and RAMClustR and method from CAMERA, so probably no
    conflicts
    - don't assign level <1 if suspect is a target? or give the choice (or make filter?)
    - spec similarity:
        - port from TPs someday
        - clarify that intensity filter happens after filtering precursor (make this optional?)
        - add proper refs
    - update filter() for new annMSMS cols
- misc
    - util to check if there are suspect results? (ie to replace inherits(...))
    - check for empty names in assertion/preparation functions
    - credits to ES
    - update version number
    - prepareSuspectList(): export?
        - mainly to allow merging of lists, perhaps make util for that instead? Would also be handy for MF databases
- expand reporting
    - eg marking which candidate corresponds to suspect and include suspect name in EICs
        - mark with different row colour and label?
    - mention suspect similarities/ranks etc for candidates (or somehow in compounds?)
    - optionally report with collapsed suspects
    - adducts from components?
- update docs & handbook
    - renamed rt/mz columns
    - new plotVenn list functionality
- tests
    - automatic InChIKey/formula calculation from InChIs/SMILES
        - already done implicitly?
    - handling empty results
    - filters
    - more?
    - new plotVenn list functionality


## docs
- improve instructions for MF and SIRIUS installation?
- improve docs for areas (only affects when features=FALSE) and average (different behavior when features=TRUE/FALSE) for as.data.table() of featureGroups
- update/check version nr mentioned in filter() for MSPeakLists
- explain xlim/ylim behavior for annotations/mols for plotSpec()


## sets
- provide methods for non-implemented functionality (eg consensus)
- find nice way to re-use docs
- filter() for features/fGroups: support ionized masses for mass filters? or just clarify it doesn't.
- handle/test empty objects
- more descriptive messages what's going on with all the avaraging of MSPeakLists
- remove sets argument for some methods (as.data.table, accessors etc)?
    - if keep, be consistent with all classes
- as.data.table() for formulas: average=T will now produce strange averaged ionized formula, for now simply remove this column.. also give a note in docs? or maybe only remove is not all adducts are equal?
- different name/generic for ionize()? makes less sense for annotation classes
- test DA algorithms
- check if more has to be cached (eg merged results from sets)
- compoundsSetMF sub-class (for settings slot)? or is access via setObjects sufficient? may need to explain anyway for other cases like intclust components
- base set class
    - don't use for fGroupsSet?
- components
    - neutralize masses? otherwise document
        - yay: consistent with other set classes
        - nay: might be a bit strange when looking for adducts etc and components are per set anyway
    - intclust: return componentsSet? if not document somewhere...
    - clearly mention that nontarget is done per set now
    - nontarget-set: plotGraph method? and make sure it's used in reportHTML()
- implement XCMS conversion functions? maybe with given set. Could just ionize() it.
- ionize() for compounds/components? if not remove formulas?
- setThreshold filter() argument, and remove argument from generators?
- handle errors when object has <=1 set
    - groupFeaturesScreening()
    - mergeScreeningSetInfos()
- suspect screening
    - handle suspect and fragment mz values somehow
    - implement TASQ?
    - consensus?
    - support recursive screening? or throw error otherwise
- fix empty MS(MS) peaklists if unavailable during merging sets
- use neutral_formula for annotated similarity calculation for formulas
- finish new minMSMSPeaks filter (apply always after averaging?)


## features
- feature optim:
    - docs
        - mention parameters default unless specified
    - keep retcor_done?
    - get rid of getXCMSSet() calls?
- suspect screening
    - rename patRoonData::targets?
    - rename groupFeaturesScreening?
- filter()
    - document which filters work on feature level (e.g. chromWidth)
    - remove zero values for maxReplicateIntRSD?
- importFeaturesXCMS/importFeaturesXCMS3/importFeatureGroupsXCMS: get rid of anaInfo arg requirement? (or make import func?)
- comparison(): support xcms3? (needs missing support for missing raw data)
- Fix: blank filter with multiple replicate groups (and maybe others?)
- Check: units of plotChord() rt/mz graphs seems off

## MSPeakLists
- isotope tagging is lost after averaging
- collapse averagedPeakLists
- test avg params
- metadata() generic?
- drop support for reAverage of [ method? doesn't seem so useful, even less so with sets


## compounds
- SIRIUS: use --auto-charge instead of manually fixing charge of fragments (or not? conflicting docs on what it does)
- test score normalization?
- timeouts for SIRIUS?
- do something about negative H explained fragments by MF?
- SusDat MF support


## formulas
- customize/document ranking column order? (only do rank for sirius?)
- getFormInfoList(): take care of consensus results like getPrecursorFormScores()

## components
- RC: check spearmans correlation
- NT: minimum size argument, combine rows for multiple rGroups?


## reporting
- add more options to reportPlots argument of reportHTML()?


## Cleanup
- Reduce non-exported class only methods


# Future

## General

- msPurity integration
- suspect screening: add MS/MS qualifiers
- fillPeaks for CAMERA (and RAMClustR?)
- support fastcluster for compounds clustering/int component clusters?
- algorithmObject() generic: for xset, xsa, rc, ...
- newProject(): fix multi line delete (when possible)
- more withr wrapping? (dev, par)
- improve default plotting for plotInt and cluster plot functions
- newProject()
    - concentration column for anaInfo
    - generate more detailed script with e.g. commented examples of subsetting, extraction etc
	- newProject(): import Bruker seq file?


## Features

- integrate OpenMS feature scoring and isotopes and PPS in general (also include filters?)
- parallel enviPick
- OpenMS MetaboliteAdductDecharger support?
- OpenMS: Support KD grouper?
- Integration of mzMine features (package pending...), MS-DIAL and KPIC2, peakonly, SIRIUS?
- suspect screening
    - tag fGroups with suspect instead of converting fGroups object (and add filter to remove non-hits)
    - automatic SMILES conversion if mz data is lacking in suspect list (and automatic name assignment of that's lacking too? might be handy for some NORMAN lists)
- topMost filter that accepts rGroups, either as AND or OR


## MSPeakLists

- DA
    - generateMSPeakListsDA: find precursor masses with larger window
    - tests
        - utils? EICs with export/vdiffr?
        - test MS peak lists deisotoping?
- metadata for Bruker peaklists?


## Formulas

- DBE calculation for SIRIUS?
- OM reporting
- as.data.table: option to average per replicate group?


## Compounds

- do something with sirius fingerprints? --> comparison?
- fix compoundViewer
- add new MF HD scorings and make sure default normalization equals that of MF web
- CFM-ID and MS-FINDER integration
- utility functions to make custom DBs for MetFrag and SIRIUS and support to use them with the latter


## components
- mass defect components
- CliqueMS
- split peak correlation and adduct etc annotation? would allow better non-target integration
- intclust
    - optionally take areas instead of intensities
    - cache results

