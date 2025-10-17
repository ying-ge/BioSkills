datacache <- new.env(hash=TRUE, parent=emptyenv())

org.Mm.eg <- function() showQCData("org.Mm.eg", datacache)
org.Mm.eg_dbconn <- function() dbconn(datacache)
org.Mm.eg_dbfile <- function() dbfile(datacache)
org.Mm.eg_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache, file=file, show.indices=show.indices)
org.Mm.eg_dbInfo <- function() dbInfo(datacache)

org.Mm.egORGANISM <- "Mus musculus"

.onLoad <- function(libname, pkgname)
{
    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "org.Mm.eg.sqlite", package=pkgname, lib.loc=libname)
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
    ann_objs <- createAnnObjs.SchemaChoice("MOUSE_DB", "org.Mm.eg", "Mouse", dbconn, datacache)
    mergeToNamespaceAndExport(ann_objs, pkgname)
    packageStartupMessage(AnnotationDbi:::annoStartupMessages("org.Mm.eg"))
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(org.Mm.eg_dbconn())
}

