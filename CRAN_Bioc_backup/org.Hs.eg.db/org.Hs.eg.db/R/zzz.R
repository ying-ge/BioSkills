datacache <- new.env(hash=TRUE, parent=emptyenv())

org.Hs.eg <- function() showQCData("org.Hs.eg", datacache)
org.Hs.eg_dbconn <- function() dbconn(datacache)
org.Hs.eg_dbfile <- function() dbfile(datacache)
org.Hs.eg_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache, file=file, show.indices=show.indices)
org.Hs.eg_dbInfo <- function() dbInfo(datacache)

org.Hs.egORGANISM <- "Homo sapiens"

.onLoad <- function(libname, pkgname)
{
    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "org.Hs.eg.sqlite", package=pkgname, lib.loc=libname)
    assign("dbfile", dbfile, envir=datacache)
    dbconn <- dbFileConnect(dbfile)
    assign("dbconn", dbconn, envir=datacache)

    ## Create the OrgDb object
    sPkgname <- sub(".db$","",pkgname)
    txdb <- loadDb(system.file("extdata", paste(sPkgname,
      ".sqlite",sep=""), package=pkgname, lib.loc=libname),
                   packageName=pkgname)    
    dbNewname <- AnnotationDbi:::dbObjectName(pkgname,"OrgDb")
    ns <- asNamespace(pkgname)
    assign(dbNewname, txdb, envir=ns)
    namespaceExport(ns, dbNewname)
        
    ## Create the AnnObj instances
    ann_objs <- createAnnObjs.SchemaChoice("HUMAN_DB", "org.Hs.eg", "Human", dbconn, datacache)
    mergeToNamespaceAndExport(ann_objs, pkgname)
    packageStartupMessage(AnnotationDbi:::annoStartupMessages("org.Hs.eg"))
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(org.Hs.eg_dbconn())
}

