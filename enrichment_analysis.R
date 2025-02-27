getPrefixRaMP = function(hmdbids = NULL, keggids = NULL, chemspiderids = NULL, CASids = NULL, chebiids = NULL, pubchemids = NULL){
  # RaMP数据库的数据准备
  return(c(na.omit(unlist(c(
    ifelse(is.null(hmdbids)|length(hmdbids)==0, NA, list(paste('hmdb', hmdbids, sep=':'))),
    ifelse(is.null(keggids)|length(keggids)==0, NA, list(paste('kegg', keggids, sep=':'))),
    ifelse(is.null(chemspiderids)|length(chemspiderids)==0, NA, list(paste('chemspider', chemspiderids, sep=':'))),
    ifelse(is.null(CASids)|length(CASids)==0, NA, list(paste('CAS', CASids, sep=':'))),
    ifelse(is.null(chebiids)|length(chebiids)==0, NA, list(paste('chebi', chebiids, sep=':'))),
    ifelse(is.null(pubchemids)|length(pubchemids)==0, NA, list(paste('pubchem', pubchemids, sep=':')))
  )))))
}

getPathEnrichPaMP = function(prefixedIDs){
  url1 = 'https://rampdb.nih.gov/api/combined-fisher-test'
  content = paste('{\n"analytes": [\n', paste('"', prefixedIDs, '"', sep = '', collapse = ',\n'), '\n]\n}', sep = ' ')
  response = httr::POST(url1, body = content, encoding = 'json')
  if (response$status_code==200){
    html = rvest::read_html(response, encoding = 'utf-8')
    htmlparse = XML::htmlTreeParse(html, useInternalNodes = T)
    txt = XML::xpathApply(htmlparse, '//body', XML::xmlValue)[[1]]
    dat = jsonlite::fromJSON(txt)
    res = dat$data$fishresults
    res = if (inherits(res, 'list') & length(res)==0) data.frame() else res
  }else{
    res = NULL
    warning(paste('Non enrichment information, code', esponse$status_code))
  }
  return(res)
}

pattEnrich  = function(pattsRes, metAnno){
  enrichRes = list()
  for (g in colnames(pattsRes)[2:length(colnames(pattsRes))]){
    ps = na.omit(unique(pattsRes[[g]]))
    # print(g)
    for (p in ps){
      onep = dplyr::filter(pattsRes, get(g) == p) 
      # print(onep)
      ids = dplyr::select(metAnno, all_of(c('COMP ID','PUBCHEM ID_nat', 'CHEMSPIDER ID_nat', 'CAS ID_nat', 'KEGG ID_nat', 'HMDB ID_nat')))
      ids = dplyr::left_join(onep, ids, by = c('metabolite' = 'COMP ID'))
      hmdbids= dplyr::filter(ids, !is.na(`HMDB ID_nat`))
      remains = dplyr::filter(ids, is.na(`HMDB ID_nat`))
      hmdbids = hmdbids$`HMDB ID_nat`
      keggids = dplyr::filter(remains, !is.na(`KEGG ID_nat`))
      remains = dplyr::filter(remains, is.na(`KEGG ID_nat`))
      keggids = keggids$`KEGG ID_nat`
      chemspiderids = dplyr::filter(remains, !is.na(`CHEMSPIDER ID_nat`))
      remains = dplyr::filter(remains, is.na(`CHEMSPIDER ID_nat`))
      chemspiderids = chemspiderids$`CHEMSPIDER ID_nat`
      casids = dplyr::filter(remains, !is.na(`CAS ID_nat`))
      remains = dplyr::filter(remains, is.na(`CAS ID_nat`))
      casids = casids$`CAS ID_nat`
      pubchemids = dplyr::filter(remains, !is.na(`PUBCHEM ID_nat`))
      remains = dplyr::filter(remains, is.na(`PUBCHEM ID_nat`))
      pubchemids = pubchemids$`PUBCHEM ID_nat`
      # print(remains)
      if (nrow(remains)>0){
        print(paste(nrow(remains), 'metabolites have no ids'))
      }
      ids = getPrefixRaMP(hmdbids = hmdbids, keggids = keggids, chemspiderids = chemspiderids,
                          CASids = casids, pubchemids = pubchemids)
      # print(ids)
      res = getPathEnrichPaMP(ids)
      # return(list(
      #   'hmdbids' = hmdbids,
      #   'keggids' = keggids,
      #   'chemspiderids' = chemspiderids,
      #   'casids' = casids
      # ))
      name = paste0(g,'+', p, collapse = '')
      enrichRes[[name]] = res
    }
  }
  return(enrichRes)
} 
