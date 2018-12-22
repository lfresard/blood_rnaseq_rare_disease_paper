#!/bin/R
#

# Functions used to generate Figures

# Functions
get_affected_status<-function(mysample, affected_status_df){
	affected_status=affected_status_df %>% filter(sample==mysample) %>% select(status)%>% pull
	return(affected_status)
}

get_institution<-function(mysample, metadata){
	institution=metadata %>% filter(sample_id==mysample) %>% select(institution)%>% pull
	return(institution)
}

get_technology<-function(mysample, metadata){
	technology=metadata %>% filter(sample_id==mysample) %>% select(variant_data)%>% pull
	return(technology)
}


get_outlier_genes<-function(mysample, threshold, outlier_df.m, pLI_thre){
	number_outliergenes=outlier_df.m %>% filter(variable==mysample)%>% filter(absZ>=threshold) %>%left_join(exac,by=c("gene"="ensgene"), all.x=T)%>% distinct %>% filter(pLI>=pLI_thre) %>%select(gene) %>% summarise(n_Elements = n_distinct(gene)) %>% pull
	return(number_outliergenes)
}

get_outlier_samples<-function(gene, threshold, outlier_df.m){
	outlier_samples=outlier.m %>% filter(gene==gene)%>% filter(abs(value)>=threshold) %>% select(variable) %>% distinct %>% pull
	return(outlier_samples)

}

get_outlier_genes_RV<-function(mysample, threshold, outlier_df.m, distance,CADD, pLI_thre){
	number_outliergenes=outlier_df.m %>% filter(variable==mysample)%>% filter(absZ>=threshold) %>% filter(distance_variant_junction<=distance)%>%
		 left_join(exac,by=c("gene_name"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=CADD, pLI>=pLI_thre) %>% 
		select(gene) %>% summarise(n_Elements = n_distinct(gene)) %>% pull
	return(number_outliergenes)
}
get_outlier_genes_RV_window<-function(mysample, threshold, outlier_df.m,CADD, pLI_thre){
	if(CADD==0){
		number_outliergenes=outlier_df.m %>% filter(variable==mysample)%>% filter(absZ>=threshold) %>%
		left_join(exac,by=c("gene_name"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=CADD| phred_cadd==".", pLI>=pLI_thre) %>% 
		select(gene) %>% summarise(n_Elements = n_distinct(gene)) %>% pull
	}
	else{
		number_outliergenes=outlier_df.m %>% filter(variable==mysample)%>% filter(absZ>=threshold) %>%
		left_join(exac,by=c("gene_name"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=CADD, pLI>=pLI_thre) %>% 
		select(gene) %>% summarise(n_Elements = n_distinct(gene)) %>% pull

	}
	return(number_outliergenes)
}

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=TRUE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=TRUE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}



grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
	plots <- list(...)
	position <- match.arg(position)
	g <- ggplotGrob(plots[[1]] + 
	theme(legend.position = position))$grobs
	legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
	lheight <- sum(legend$height)
	lwidth <- sum(legend$width)
	gl <- lapply(plots, function(x) x +
	theme(legend.position = "none"))
	gl <- c(gl, ncol = ncol, nrow = nrow)

	combined <- switch(position,
	                   "bottom" = arrangeGrob(do.call(arrangeGrob, gl), 
	                   legend,ncol = 1,
					heights = unit.c(unit(1, "npc") - lheight, lheight)),
					"right" = arrangeGrob(do.call(arrangeGrob, gl),
				  legend, ncol = 2,
					widths = unit.c(unit(1, "npc") - lwidth, lwidth)))

	grid.newpage()
	grid.draw(combined)

	# return gtable invisibly
	invisible(combined)
}

process_variant=function(variant_freq){
	variant_freq=as.character(variant_freq)
	if (grepl(",", variant_freq)) {
		variant_freq=max(as.numeric(unlist(strsplit(variant_freq,","))))
		return(as.numeric(variant_freq))
	}
	else{ return(as.numeric(variant_freq))}
}


get_outlier_samples_signed<-function(mygene, threshold, outlier_df.m, method, affected_status_df){
	if (method == "<"){
		outliers=outlier_df.m%>% filter(gene==mygene)%>% filter(value<=threshold) %>% select(variable) %>% distinct
		outliers$variable=factor(outliers$variable, levels=levels(affected_status_df$sample))
	}else if (method == ">"){
		outliers=outlier_df.m%>% filter(gene==mygene)%>% filter(value>=threshold) %>% select(variable) %>% distinct
		outliers$variable=factor(outliers$variable, levels=levels(affected_status_df$sample))
	}else print("Error: choose '<' or '>'")
			
	if(nrow(outliers)==0){
		outliers_counts=data.frame(Case=0, Control=0, GeneID=mygene)
		return(outliers_counts)
	}else{
		outliers$affected_status=sapply(outliers$variable, get_affected_status,affected_status_df=affected_status_df)
		outlier_counts=outliers%>% group_by(affected_status) %>% summarise(n=n())  %>% complete(affected_status, fill=list(n=0)) %>% as.data.frame %>% spread(affected_status, n)
		outlier_counts$GeneID =mygene
		return(outlier_counts)
}}


get_outlier_samples_signed_abs<-function(mygene, threshold, outlier_df.m, method, affected_status_df){
	if (method == "<"){
		outliers=outlier_df.m%>% filter(gene==mygene)%>% filter(value<=threshold) %>% select(variable) %>% distinct
		outliers$variable=factor(outliers$variable, levels=levels(affected_status_df$sample))
	}else if (method == ">"){
		outliers=outlier_df.m%>% filter(gene==mygene)%>% filter(abs(value)>=threshold) %>% select(variable) %>% distinct
		outliers$variable=factor(outliers$variable, levels=levels(affected_status_df$sample))
	}else print("Error: choose '<' or '>'")
			
	if(nrow(outliers)==0){
		outliers_counts=data.frame(Case=0, Control=0, GeneID=mygene)
		return(outliers_counts)
	}else{
		outliers$affected_status=sapply(outliers$variable, get_affected_status,affected_status_df=affected_status_df)
		outlier_counts=outliers%>% group_by(affected_status) %>% summarise(n=n())  %>% complete(affected_status, fill=list(n=0)) %>% as.data.frame %>% spread(affected_status, n)
		outlier_counts$GeneID =mygene
		return(outlier_counts)
}}
is_annotated<-function(junc_ref_vect, junc_test){
	annot=ifelse(junc_test %in% junc_ref_vect, "A", "NA")
	return(annot)
}

get_corrected_data<-function(affected_status_df, ratio_df){
	#filter for samples of interest
	samples=affected_status_df$sample
	ratio_df=ratio_df[,samples]
	#pca correction
	pca=prcomp(t(ratio_df))
	sum_pca=t(summary(pca)$importance)
	colnames(sum_pca)=c("sd", "prop_var", "cumul_var")
	pcs_number=as.data.frame(sum_pca) %>% filter(cumul_var <=0.95) %>% nrow
	pcs_selec=pca$x[,1:pcs_number]
	mod = model.matrix(~1, data=as.data.frame(t(ratio_df)))
	modsv <- cbind(mod, pcs_selec)
	fitsv <- lm.fit(modsv, t(ratio_df))


	junction_zscore_matrix=as.matrix(fitsv$residuals)
	
	## Scale and center
	junction_zscore_scale =scale(junction_zscore_matrix, center=TRUE, scale=TRUE)
	junction_zscore_scale.t=t(junction_zscore_scale)
	junction_zscore_scale.t.df=as.data.frame(junction_zscore_scale.t)
	
	junction_zscore_scale.t.df$annotation_status=sapply(junctions,is_annotated, junc_ref_vect=annot_junctions)
	junction_zscore_scale.t.df$chr=unlist(strsplit(junctions, "_"))[c(TRUE,FALSE,FALSE,FALSE)]
	junction_zscore_scale.t.df$junc_start=unlist(strsplit(junctions, "_"))[c(FALSE,TRUE,FALSE,FALSE)]
	junction_zscore_scale.t.df$junc_end=unlist(strsplit(junctions, "_"))[c(FALSE,FALSE,TRUE,FALSE)]
	junction_zscore_scale.t.df$gene=unlist(strsplit(junctions, "_"))[c(FALSE,FALSE,FALSE, TRUE)]

	
	junction_zscore_scale.t.df.m=melt(junction_zscore_scale.t.df, id.var=c("chr", "junc_start", "junc_end","gene", "annotation_status"))
	junction_zscore_scale.t.df.m$absZ=abs(junction_zscore_scale.t.df.m$value)
	junction_zscore_scale.t.df.m=junction_zscore_scale.t.df.m %>% left_join(affected_status_df, by=c("variable"="sample")) %>% select(chr ,junc_start  ,junc_end, gene ,annotation_status ,variable, value, absZ, status)

	return(junction_zscore_scale.t.df.m)
}


process_variant=function(variant_freq){
	variant_freq=as.character(variant_freq)
	if (grepl(",", variant_freq)) {
		variant_freq=max(as.numeric(unlist(strsplit(variant_freq,","))))
		return(variant_freq)
	}
	else{ return(as.numeric(variant_freq))}
}

get_number_junction_RV=function(sample,sample_dir, MAF, mydistance, sgl){
	file_name=paste0(sample_dir, sample,"_junction_RV.txt")
	temp=as.data.frame(read.table(file_name))
	colnames(temp)=c("chr_var", "start_var", "end_var", "ref", "alt", "gnomAD_AF","raw_cadd", "phred_cadd", "chr_junc", "start_junc","end_junc", "gene_junc","chr_gene", "start_gene","end_gene", "gene_name", "distance")
	if(sgl=="yes"){
		temp= temp %>% filter(distance<=mydistance) %>%mutate(gnomAD_AF = replace(gnomAD_AF, gnomAD_AF==".", 0)) 
		}else{temp= temp %>% filter(distance>=mydistance) %>%mutate(gnomAD_AF = replace(gnomAD_AF, gnomAD_AF==".", NA)) }
	temp$gnomAD_AF_filt=sapply(temp$gnomAD_AF,process_variant)
	RV_number=temp %>% filter(gnomAD_AF_filt <= MAF)%>% nrow 
	return(RV_number)
}
get_number_junction_RV_window=function(sample,sample_dir, MAF, sgl){
	file_name=paste0(sample_dir, sample,"_junction_RV_window.txt")
	temp=as.data.frame(read.table(file_name))
	colnames(temp)=c("chr_var", "start_var", "end_var", "ref", "alt", "gnomAD_AF","raw_cadd", "phred_cadd", "chr_junc", "start_junc","end_junc", "gene_junc","chr_gene", "start_gene","end_gene", "gene_name")
	if(sgl=="yes"){
		temp= temp%>%mutate(gnomAD_AF = replace(gnomAD_AF, gnomAD_AF==".", 0)) 
		}else{temp= temp %>%mutate(gnomAD_AF = replace(gnomAD_AF, gnomAD_AF==".", NA)) }
	temp$gnomAD_AF_filt=sapply(temp$gnomAD_AF,process_variant)
	RV_number=temp %>% filter(gnomAD_AF_filt <= MAF)%>% nrow 
	return(RV_number)
}

get_number_junction_RVgene=function(sample,sample_dir, MAF, sgl, CADD){
	file_name=paste0(sample_dir, sample,"_junction_RV_window.txt")
	temp=as.data.frame(read.table(file_name))
	colnames(temp)=c("chr_var", "start_var", "end_var", "ref", "alt", "gnomAD_AF","raw_cadd", "phred_cadd", "chr_junc", "start_junc","end_junc", "gene_junc","chr_gene", "start_gene","end_gene", "gene_name")
	if(sgl=="yes"){
		temp= temp %>% filter(as.numeric(as.character(phred_cadd))>=CADD) %>%mutate(gnomAD_AF = replace(gnomAD_AF, gnomAD_AF==".", 0)) 
		}else{temp= temp %>% filter(as.numeric(as.character(phred_cadd))>=CADD) %>%mutate(gnomAD_AF = replace(gnomAD_AF, gnomAD_AF==".", NA)) }
	temp$gnomAD_AF_filt=sapply(temp$gnomAD_AF,process_variant)
	gene_RV_number=temp %>% filter(gnomAD_AF_filt <= MAF)%>%  select(gene_name) %>% unique %>%nrow  
	return(gene_RV_number)
}
get_number_junction_RVgene_causalgene=function(sample,sample_dir, MAF, sgl, CADD, causal_gene_name){
	file_name=paste0(sample_dir, sample,"_junction_RV_window.txt")
	temp=as.data.frame(read.table(file_name))
	colnames(temp)=c("chr_var", "start_var", "end_var", "ref", "alt", "gnomAD_AF","raw_cadd", "phred_cadd", "chr_junc", "start_junc","end_junc", "gene_junc","chr_gene", "start_gene","end_gene", "gene_name")
	if(sgl=="yes"){
		temp= temp %>% filter(as.numeric(as.character(phred_cadd))>=CADD) %>%mutate(gnomAD_AF = replace(gnomAD_AF, gnomAD_AF==".", 0)) 
		}else{temp= temp %>% filter(as.numeric(as.character(phred_cadd))>=CADD) %>%mutate(gnomAD_AF = replace(gnomAD_AF, gnomAD_AF==".", NA)) }
	temp$gnomAD_AF_filt=sapply(temp$gnomAD_AF,process_variant)
	gene_RV=temp %>% filter(gnomAD_AF_filt <= MAF)%>%  select(gene_name) %>% mutate(gene_name=str_extract(gene_name, "ENSG[0-9]+"))%>% unique  
	if(causal_gene_name %in% gene_RV$gene_name){
		return(1)
	} else{return(0)}
}

get_number_junction_RVgeneHPO=function(sample,sample_dir, MAF, sgl, CADD){
	file_name=paste0(sample_dir, sample,"_junction_RV_window.txt")
	temp=as.data.frame(read.table(file_name))
	colnames(temp)=c("chr_var", "start_var", "end_var", "ref", "alt", "gnomAD_AF","raw_cadd", "phred_cadd", "chr_junc", "start_junc","end_junc", "gene_junc","chr_gene", "start_gene","end_gene", "gene_name")
	if(sgl=="yes"){
		temp= temp %>% filter(as.numeric(as.character(phred_cadd))>=CADD) %>%mutate(gnomAD_AF = replace(gnomAD_AF, gnomAD_AF==".", 0)) %>% mutate(gene_name=str_extract(gene_name, "ENSG[0-9]+"))%>%left_join(grch37, by=c("gene_name"="ensgene")) %>%select(gene_name, symbol, entrez, phred_cadd,gnomAD_AF)%>%  left_join(gene_to_pheno, by=c("entrez", "symbol")) %>% distinct %>% left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>% mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPO_extended[[sample]], 1, 0))  %>% select(gene_name, symbol,gnomAD_AF, phred_cadd, hpo_in_input) %>% distinct %>% filter(hpo_in_input>=1) %>% distinct
		}else{temp= temp %>% filter(as.numeric(as.character(phred_cadd))>=CADD) %>%mutate(gnomAD_AF = replace(gnomAD_AF, gnomAD_AF==".", NA))%>% mutate(gene_name=str_extract(gene_name, "ENSG[0-9]+"))%>%left_join(grch37, by=c("gene_name"="ensgene")) %>%select(gene_name, symbol, entrez, phred_cadd,gnomAD_AF)%>%  left_join(gene_to_pheno, by=c("entrez", "symbol")) %>% distinct %>% left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>% mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPO_extended[[sample]], 1, 0))  %>% select(gene_name, symbol,gnomAD_AF, phred_cadd, hpo_in_input) %>% distinct %>% filter(hpo_in_input>=1) %>% distinct }
	temp$gnomAD_AF_filt=sapply(temp$gnomAD_AF,process_variant)
	gene_RV_number=temp %>% filter(gnomAD_AF_filt <= MAF)%>%  select(gene_name) %>% unique %>%nrow  
	return(gene_RV_number)
}

get_number_junction_RVgeneHPOcausalgene=function(sample,sample_dir, MAF, sgl, CADD, causal_gene_name){
	file_name=paste0(sample_dir, sample,"_junction_RV_window.txt")
	temp=as.data.frame(read.table(file_name))
	colnames(temp)=c("chr_var", "start_var", "end_var", "ref", "alt", "gnomAD_AF","raw_cadd", "phred_cadd", "chr_junc", "start_junc","end_junc", "gene_junc","chr_gene", "start_gene","end_gene", "gene_name")
	if(sgl=="yes"){
		temp= temp %>% filter(as.numeric(as.character(phred_cadd))>=CADD) %>%mutate(gnomAD_AF = replace(gnomAD_AF, gnomAD_AF==".", 0)) %>% mutate(gene_name=str_extract(gene_name, "ENSG[0-9]+"))%>%left_join(grch37, by=c("gene_name"="ensgene")) %>%select(gene_name, symbol, entrez, phred_cadd,gnomAD_AF)%>%  left_join(gene_to_pheno, by=c("entrez", "symbol")) %>% distinct %>% left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>% mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPO_extended[[sample]], 1, 0))  %>% select(gene_name, symbol,gnomAD_AF, phred_cadd, hpo_in_input) %>% distinct %>% filter(hpo_in_input>=1) %>% distinct
		}else{temp= temp %>% filter(as.numeric(as.character(phred_cadd))>=CADD) %>%mutate(gnomAD_AF = replace(gnomAD_AF, gnomAD_AF==".", NA))%>% mutate(gene_name=str_extract(gene_name, "ENSG[0-9]+"))%>%left_join(grch37, by=c("gene_name"="ensgene")) %>%select(gene_name, symbol, entrez, phred_cadd,gnomAD_AF)%>%  left_join(gene_to_pheno, by=c("entrez", "symbol")) %>% distinct %>% left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>% mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPO_extended[[sample]], 1, 0))  %>% select(gene_name, symbol,gnomAD_AF, phred_cadd, hpo_in_input) %>% distinct %>% filter(hpo_in_input>=1) %>% distinct }
	temp$gnomAD_AF_filt=sapply(temp$gnomAD_AF,process_variant)
	gene_RV=temp %>% filter(gnomAD_AF_filt <= MAF)%>%  select(gene_name) %>% unique 
	if(causal_gene_name %in% gene_RV$gene_name){
		return(1)
	} else{return(0)}
}


get_number_RVgene_Hpo=function(sample,sample_dir, MAF,sgl, CADD){
	file_name=paste0(sample_dir, sample,"_homogenized_gnomad_cadd_RV_withgene.bed.gz")
	temp=as.data.frame(read.table(gzfile(file_name)))
	colnames(temp)=c("chr_var", "start_var", "end_var", "ref", "alt", "gnomAD_AF","raw_cadd", "phred_cadd", "chr_gene", "start_gene","end_gene", "gene_name")
	if(sgl=="yes"){
		temp= temp %>% filter(as.numeric(as.character(phred_cadd))>=CADD)%>% mutate(gnomAD_AF = replace(gnomAD_AF, gnomAD_AF==".", 0))%>% mutate(gene_name=str_extract(gene_name, "ENSG[0-9]+"))%>%left_join(grch37, by=c("gene_name"="ensgene")) %>%select(gene_name, symbol, entrez, phred_cadd,gnomAD_AF)%>%  left_join(gene_to_pheno, by=c("entrez", "symbol")) %>% distinct %>% left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>% mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPO_extended[[sample]], 1, 0)) %>% select(gene_name, symbol,gnomAD_AF, phred_cadd, hpo_in_input) %>% distinct %>% filter(hpo_in_input>=1) %>% distinct
		}else{temp= temp%>% filter(as.numeric(as.character(phred_cadd))>=CADD) %>% mutate(gnomAD_AF = replace(gnomAD_AF, gnomAD_AF==".", NA))%>% mutate(gene_name=str_extract(gene_name, "ENSG[0-9]+"))%>%left_join(grch37, by=c("gene_name"="ensgene")) %>%select(gene_name, symbol, entrez, phred_cadd,gnomAD_AF)%>% left_join(gene_to_pheno, by=c("entrez", "symbol")) %>% distinct %>% left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>% mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPO_extended[[sample]], 1, 0)) %>% select(gene_name, symbol,gnomAD_AF, phred_cadd, hpo_in_input) %>% distinct %>% filter(hpo_in_input>=1) %>% distinct }
	temp$gnomAD_AF_filt=sapply(temp$gnomAD_AF,process_variant)
	gene_RV_number=temp %>% filter(gnomAD_AF_filt <= MAF)%>% select(gene_name) %>% unique %>%nrow 
	return(gene_RV_number)
}

get_number_RVgene_Hpo_causalgene=function(sample,sample_dir, MAF,sgl, CADD, causal_gene_name){
	file_name=paste0(sample_dir, sample,"_homogenized_gnomad_cadd_RV_withgene.bed.gz")
	temp=as.data.frame(read.table(gzfile(file_name)))
	colnames(temp)=c("chr_var", "start_var", "end_var", "ref", "alt", "gnomAD_AF","raw_cadd", "phred_cadd", "chr_gene", "start_gene","end_gene", "gene_name")
	if(sgl=="yes"){
		temp= temp %>% filter(as.numeric(as.character(phred_cadd))>=CADD)%>% mutate(gnomAD_AF = replace(gnomAD_AF, gnomAD_AF==".", 0))%>% mutate(gene_name=str_extract(gene_name, "ENSG[0-9]+"))%>%left_join(grch37, by=c("gene_name"="ensgene")) %>%select(gene_name, symbol, entrez, phred_cadd,gnomAD_AF)%>%  left_join(gene_to_pheno, by=c("entrez", "symbol")) %>% distinct %>% left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>% mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPO_extended[[sample]], 1, 0)) %>% select(gene_name, symbol,gnomAD_AF, phred_cadd, hpo_in_input) %>% distinct %>% filter(hpo_in_input>=1) %>% distinct
		}else{temp= temp%>% filter(as.numeric(as.character(phred_cadd))>=CADD) %>% mutate(gnomAD_AF = replace(gnomAD_AF, gnomAD_AF==".", NA))%>% mutate(gene_name=str_extract(gene_name, "ENSG[0-9]+"))%>%left_join(grch37, by=c("gene_name"="ensgene")) %>%select(gene_name, symbol, entrez, phred_cadd,gnomAD_AF)%>% left_join(gene_to_pheno, by=c("entrez", "symbol")) %>% distinct %>% left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>% mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPO_extended[[sample]], 1, 0)) %>% select(gene_name, symbol,gnomAD_AF, phred_cadd, hpo_in_input) %>% distinct %>% filter(hpo_in_input>=1) %>% distinct }
	temp$gnomAD_AF_filt=sapply(temp$gnomAD_AF,process_variant)
	gene_RV=temp %>% filter(gnomAD_AF_filt <= MAF)%>% select(gene_name) %>% unique
	if(causal_gene_name %in% gene_RV$gene_name){
		return(1)
	} else{return(0)}
}


get_number_RV=function(sample,sample_dir, MAF,sgl){
	file_name=paste0(sample_dir, sample,"_homogenized_gnomad_cadd_RV_withgene.bed.gz")
	temp=as.data.frame(read.table(gzfile(file_name)))
	colnames(temp)=c("chr_var", "start_var", "end_var", "ref", "alt", "gnomAD_AF","raw_cadd", "phred_cadd", "chr_gene", "start_gene","end_gene", "gene_name")
	if(sgl=="yes"){
		temp= temp %>% mutate(gnomAD_AF = replace(gnomAD_AF, gnomAD_AF==".", 0))
		}else{temp= temp %>% mutate(gnomAD_AF = replace(gnomAD_AF, gnomAD_AF==".", NA)) }
	temp$gnomAD_AF_filt=sapply(temp$gnomAD_AF,process_variant)
	RV_number=temp %>% filter(gnomAD_AF_filt <= MAF)%>% nrow 
	return(RV_number)
}


# average replaced by median to see if eliminates saples with too many variants
get_number_RV_pergene=function(sample,sample_dir, MAF,sgl, gene_list){
	file_name=paste0(sample_dir, sample,"_homogenized_gnomad_cadd_RV_withgene.bed.gz")
	temp=as.data.frame(read.table(gzfile(file_name)))
	colnames(temp)=c("chr_var", "start_var", "end_var", "ref", "alt", "gnomAD_AF","raw_cadd", "phred_cadd", "chr_gene", "start_gene","end_gene", "gene_name")
	if(sgl=="yes"){
		temp= temp %>% mutate(gnomAD_AF = replace(gnomAD_AF, gnomAD_AF==".", 0))
		}else{temp= temp %>% mutate(gnomAD_AF = replace(gnomAD_AF, gnomAD_AF==".", NA)) }
	temp$gnomAD_AF_filt=sapply(temp$gnomAD_AF,process_variant)
	gene_vector=gene_list[[sample]][,1]
	RV_number= temp %>% mutate(gene_name=str_extract(gene_name, "ENSG[0-9]+")) %>% filter(gene_name %in% gene_vector) %>% filter(gnomAD_AF_filt <= MAF)%>% filter(as.numeric(as.character(phred_cadd))>=10) %>% group_by(gene_name) %>% summarize(n_var=n()) %>% ungroup %>% summarize(median=median(n_var, na.rm=T)) %>% pull
	return(RV_number)
}


get_number_RVgene=function(sample,sample_dir, MAF,sgl, CADD){
	file_name=paste0(sample_dir, sample,"_homogenized_gnomad_cadd_RV_withgene.bed.gz")
	temp=as.data.frame(read.table(gzfile(file_name)))
	colnames(temp)=c("chr_var", "start_var", "end_var", "ref", "alt", "gnomAD_AF","raw_cadd", "phred_cadd", "chr_gene", "start_gene","end_gene", "gene_name")
	if(sgl=="yes"){
		temp= temp %>% filter(as.numeric(as.character(phred_cadd))>=CADD)%>% mutate(gnomAD_AF = replace(gnomAD_AF, gnomAD_AF==".", 0))
		}else{temp= temp%>% filter(as.numeric(as.character(phred_cadd))>=CADD) %>% mutate(gnomAD_AF = replace(gnomAD_AF, gnomAD_AF==".", NA)) }
	temp$gnomAD_AF_filt=sapply(temp$gnomAD_AF,process_variant)
	gene_RV_number=temp %>% filter(gnomAD_AF_filt <= MAF)%>% select(gene_name) %>% filter (! duplicated(gene_name)) %>%nrow 
	return(gene_RV_number)
}

get_number_RVcausalgene=function(sample,sample_dir, MAF,sgl, CADD, causal_gene_name){
	file_name=paste0(sample_dir, sample,"_homogenized_gnomad_cadd_RV_withgene.bed.gz")
	temp=as.data.frame(read.table(gzfile(file_name)))
	colnames(temp)=c("chr_var", "start_var", "end_var", "ref", "alt", "gnomAD_AF","raw_cadd", "phred_cadd", "chr_gene", "start_gene","end_gene", "gene_name")
	if(sgl=="yes"){
		temp= temp %>% filter(as.numeric(as.character(phred_cadd))>=CADD)%>% mutate(gnomAD_AF = replace(gnomAD_AF, gnomAD_AF==".", 0))
		}else{temp= temp%>% filter(as.numeric(as.character(phred_cadd))>=CADD) %>% mutate(gnomAD_AF = replace(gnomAD_AF, gnomAD_AF==".", NA)) }
	temp$gnomAD_AF_filt=sapply(temp$gnomAD_AF,process_variant)
	gene_RV=temp %>% filter(gnomAD_AF_filt <= MAF)%>% select(gene_name) %>% mutate(gene_name=unlist(strsplit(as.character(gene_name),"[.]"))[c(TRUE,FALSE)])%>%unique  
	if(causal_gene_name %in% gene_RV$gene_name){
		return(1)
	} else{return(0)}
}


get_number_RV_outlier=function(sample,outlier_file,zscore_threshold, MAF, mydistance, sgl){
	#outier_file= outier_file %>% mutate(gnomAD_AF_withsgl = replace(gnomAD_AF, gnomAD_AF==".", 0))
	if(sgl=="yes"){
		out_sample_number=outlier_file %>% filter(variable==sample, absZ>=zscore_threshold, gnomAD_AF_filt_withsgl<= MAF, distance_variant_junction<=mydistance)%>% select(chr_var,start_var,end_var) %>% distinct%>% nrow
	}else{ 	out_sample_number=outlier_file %>% filter(variable==sample, absZ>=zscore_threshold, gnomAD_AF_filt_nosgl<= MAF, distance_variant_junction<=mydistance)%>% select(chr_var,start_var,end_var) %>% distinct%>% nrow
}
	return(out_sample_number)
}



get_number_RV_outlier_window=function(sample,outlier_file,zscore_threshold, MAF, sgl){
	#outier_file= outier_file %>% mutate(gnomAD_AF_withsgl = replace(gnomAD_AF, gnomAD_AF==".", 0))
	if(sgl=="yes"){
		out_sample_number=outlier_file %>% filter(variable==sample, absZ>=zscore_threshold, gnomAD_AF_filt_withsgl<= MAF)%>% select(chr_var,start_var,end_var) %>% distinct%>% nrow
	}else{ 	out_sample_number=outlier_file %>% filter(variable==sample, absZ>=zscore_threshold, gnomAD_AF_filt_nosgl<= MAF)%>% select(chr_var,start_var,end_var) %>% distinct%>% nrow
}
	return(out_sample_number)
}

get_genelist_prop<-function(mysample,threshold, outlier_df.m,  gene_list){
	genes=outlier_df.m %>% filter(variable==mysample)%>%filter(absZ>=threshold)  %>%select(gene) 
	if(nrow(genes)!=0){
		genes=genes%>% mutate(gene=unlist(strsplit(as.character(gene),"[.]"))[c(TRUE,FALSE)]) %>% distinct %>% pull
		prop_gene_list=length(which(genes %in% gene_list))/length(genes)
		return(prop_gene_list)}
	else{return(0)
	}
}

get_genelist_prop_RV<-function(mysample,threshold, outlier_df.m, distance,cadd, gene_list){
	genes=outlier_df.m %>% filter(variable==mysample)%>%filter(absZ>=threshold) %>% filter(distance_variant_junction<=distance) %>% filter(as.numeric(as.character(phred_cadd))>=cadd)%>%select(gene) 
	if(nrow(genes)!=0){
		genes=genes%>% mutate(gene=unlist(strsplit(as.character(gene),"[.]"))[c(TRUE,FALSE)]) %>% distinct %>% pull
		prop_gene_list=length(which(genes %in% gene_list))/length(genes)
		return(prop_gene_list)
	}
	else{return(0)}
}


read_in_disease_lists<-function( disease_file,path_to_file){
	temp=read.table(paste0(path_to_file, disease_file), header=F)[,1]
	return(temp)
}


check_gene_names<-function(gene_list){
	if(all(startsWith(as.character(gene_list), "ENSG")) == TRUE){
		return(gene_list)
	}else {

		gene_list2=as.data.frame(gene_list) 
		names(gene_list2)=c('symbol')
		gene_list2=gene_list2 %>% merge(grch37,  by='symbol') %>% select(ensgene) %>% distinct
		gene_list2=gene_list2[,1]
		return(gene_list2)
	}
}

get_pct_disease_blood<-function(disease_gene_list, pct_df){
	subset=pct_df[pct_df$GENE %in% disease_gene_list,]
	res=colSums(subset[,2:ncol(subset)])
	res=res/nrow(subset)*100
	return(res)
}


make_bins=function(df, myrange){
	transform(df, bin = cut(avg_tpm, myrange, right=F))
}

make_bins_median=function(df, myrange){
	transform(df, bin = cut(median_tpm, myrange, right=F))
}


get_pct_disease_blood_bin<-function(disease_gene_list, pct_df){
	subset=pct_df[pct_df$GENE %in% disease_gene_list,]
	subset=subset %>% group_by(bin)%>% mutate( pct_genes=n()/length(disease_gene_list)*100)%>% select(bin,  pct_genes)%>% distinct
	
	return(as.data.frame(subset))
}


get_index<-function(bin_value, bin_df){
	index=bin_df$INDEX[bin_df$bin==bin_value]
	return(index)
}



proportion.ratios = function(zscore_table, thresholds) {
    # make empty data frame and fill it from rbinds (not very many, so it's fine)
    results = data.frame(ESTIM=numeric(), CI.LOW=numeric(), CI.HIGH=numeric(), THRESH = numeric(), COUNT = numeric(),
        stringsAsFactors = F)
    for (thresh in thresholds) {
        props = proportion.ratios.helper(zscore_table, thresh)
        if (!is.null(props)) {
            results = rbind(results, props)
        }
    }
    return(results)
}

proportion.ratios.helper = function(counts, thresh) {
	test= counts %>% mutate(outlier=ifelse(abs(value)>=thresh, 1,0)) %>% select(gene, variable, disease_gene, outlier) %>% filter(!is.na(outlier)) %>% group_by(gene,variable) %>% summarize(disease_gene=max(disease_gene), outlier=max(outlier))%>% ungroup  %>% select(gene,disease_gene,outlier)
	test.out=test %>% group_by(gene, disease_gene) %>% summarize(outlier=max(outlier, na.rm=T)) %>% ungroup %>% select(disease_gene,outlier)
	test.glm=glm( test.out$disease_gene~test.out$outlier, family=binomial)
	mod_summary =summary(test.glm)
	p.value=mod_summary[[12]][2, 4] # pvalue
	#p.value=coef(mod_summary)[2,4] # pvalue
	log_odds=mod_summary[[12]][2, 1] # coefficient estimate
	#log_odds=coef(mod_summary)[2, 1] # coefficient estimate
	min.ci = mod_summary[[12]][2, 1] - (1.96*mod_summary[[12]][2, 2]) #
	#min.ci = coef(mod_summary)[2, 1] - (1.96*coef(mod_summary)[2, 2]) #
	max.ci = mod_summary[[12]][2, 1] + (1.96*mod_summary[[12]][2, 2]) #
	#max.ci = coef(mod_summary)[2, 1] + (1.96*coef(mod_summary)[2, 2]) #

	#summary.counts = as.data.frame(table(test))
    #if (nrow(summary.counts) != 4 | min(summary.counts$Freq) == 0) {
    #    cat("Warning: Skipping this threshold because zeros in contingency table.\n")
    #    return(NULL)
    #}
	#out.disease = summary.counts$Freq[summary.counts$disease_gene == 1 & summary.counts$outlier == 1]
	#nonout.disease = summary.counts$Freq[summary.counts$disease_gene == 1 & summary.counts$outlier == 0]
	#out.total = sum(summary.counts$Freq[summary.counts$outlier == 1])
	#nonout.total = sum(summary.counts$Freq[summary.counts$outlier == 0])
	#estimate = (out.disease/out.total)/(nonout.disease/nonout.total)
	## get bounds of confidence interval on the log of the proportion then exponentiate
	#log.se = sqrt(1/out.disease - 1/out.total + 1/nonout.disease - 1/nonout.total)
	#max.ci = estimate * exp(1.96*log.se)
	#min.ci = estimate * exp(-1.96*log.se)
	dfrow = list(LOG_ODDS=log_odds, CI.LOW=min.ci, CI.HIGH=max.ci, THRESH=thresh, COUNT = sum(test.out$outlier), PVALUE=p.value)
    return (dfrow)
}

proportion.ratios.exp.signed = function(zscore_table, thresholds) {
    # make empty data frame and fill it from rbinds (not very many, so it's fine)
    results = data.frame(ESTIM=numeric(), CI.LOW=numeric(), CI.HIGH=numeric(), THRESH = numeric(), COUNT = numeric(),
        stringsAsFactors = F)
    for (thresh in thresholds) {
        props = proportion.ratios.exp.signed.helper (zscore_table, thresh)
        if (!is.null(props)) {
            results = rbind(results, props)
        }
    }
    return(results)
}
proportion.ratios.exp = function(zscore_table, thresholds) {
    # make empty data frame and fill it from rbinds (not very many, so it's fine)
    results = data.frame(ESTIM=numeric(), CI.LOW=numeric(), CI.HIGH=numeric(), THRESH = numeric(), COUNT = numeric(),
        stringsAsFactors = F)
    for (thresh in thresholds) {
        props = proportion.ratios.exp.helper (zscore_table, thresh)
        if (!is.null(props)) {
            results = rbind(results, props)
        }
    }
    return(results)
}


proportion.ratios.exp.signed.helper = function(counts, thresh) {
	test= counts %>% mutate(outlier=ifelse(zscore<=-thresh, 1,0)) %>% select(gene, disease_gene, outlier) 
	test.out=test %>% group_by(gene, disease_gene) %>% summarize(outlier=max(outlier)) %>% ungroup %>% select(disease_gene,outlier)
	test.glm=glm(test.out$disesease_gene ~test.out$outlier, family=binomial)
	mod_summary =summary(test.glm)
	p.value=mod_summary[[12]][2, 4] # pvalue
	log_odds=mod_summary[[12]][2, 1] # coefficient estimate
	min.ci = mod_summary[[12]][2, 1] - (1.96*mod_summary[[12]][2, 2]) #
	max.ci = mod_summary[[12]][2, 1] + (1.96*mod_summary[[12]][2, 2]) #

	dfrow = list(LOG_ODDS=log_odds, CI.LOW=min.ci, CI.HIGH=max.ci, THRESH=thresh, COUNT = sum(test.out$outlier), PVALUE=p.value)
    return (dfrow)
}


proportion.ratios.exp.helper = function(counts, thresh) {
	test= counts %>% mutate(outlier=ifelse(abs(zscore)>=thresh, 1,0)) %>% select(gene, disease_gene, outlier) 
	test.out=test %>% group_by(gene, disease_gene) %>% summarize(outlier=max(outlier)) %>% ungroup %>% select(disease_gene,outlier)
	test.glm=glm(test.out$disease_gene ~test.out$outlier, family=binomial)
	mod_summary =summary(test.glm)
	p.value=mod_summary[[12]][2, 4] # pvalue
	log_odds=mod_summary[[12]][2, 1] # coefficient estimate
	min.ci = mod_summary[[12]][2, 1] - (1.96*mod_summary[[12]][2, 2]) #
	max.ci = mod_summary[[12]][2, 1] + (1.96*mod_summary[[12]][2, 2]) #
	


	#summary.counts = as.data.frame(table(test))
    #if (nrow(summary.counts) != 4 | min(summary.counts$Freq) == 0) {
    #    cat("Warning: Skipping this threshold because zeros in contingency table.\n")
    #    return(NULL)
    #}
	#out.disease = summary.counts$Freq[summary.counts$disease_gene == 1 & summary.counts$outlier == 1]
	#nonout.disease = summary.counts$Freq[summary.counts$disease_gene == 1 & summary.counts$outlier == 0]
	#out.total = sum(summary.counts$Freq[summary.counts$outlier == 1])
	#nonout.total = sum(summary.counts$Freq[summary.counts$outlier == 0])
	#estimate = (out.disease/out.total)/(nonout.disease/nonout.total)
	## get bounds of confidence interval on the log of the proportion then exponentiate
	#log.se = sqrt(1/out.disease - 1/out.total + 1/nonout.disease - 1/nonout.total)
	#max.ci = estimate * exp(1.96*log.se)
	#min.ci = estimate * exp(-1.96*log.se)
	dfrow = list(LOG_ODDS=log_odds, CI.LOW=min.ci, CI.HIGH=max.ci, THRESH=thresh, COUNT = sum(test.out$outlier), PVALUE=p.value)
    return (dfrow)
}


proportion.ratios.combined = function(spli_zscores, exp_zscores, thresholds) {
    # make empty data frame and fill it from rbinds (not very many, so it's fine)
    results = data.frame(LOG_ODDS=numeric(), CI.LOW=numeric(), CI.HIGH=numeric(), THRESH = numeric(), COUNT = numeric(),PVALUE=numeric(),
        stringsAsFactors = F)
    for (thresh in thresholds) {
        props = proportion.ratios.combined.helper(spli_zscores,exp_zscores, thresh)
        if (!is.null(props)) {
            results = rbind(results, props)
        }
    }
    return(results)
}

proportion.ratios.combined.helper = function(spli.counts,exp.counts, thresh) {
	test.spli= spli.counts %>% mutate(spli.outlier=ifelse(abs(value)>=thresh, 1,0)) %>% select(gene, variable, disease_gene, spli.outlier) %>% group_by(gene,variable) %>% summarize(disease_gene=max(disease_gene), spli.outlier=max(spli.outlier))%>% ungroup  %>% select(gene, variable,disease_gene, spli.outlier)%>% rename(sample_id=variable)
	test.exp=exp.counts%>% mutate(exp.outlier=ifelse(abs(zscore)>=thresh, 1,0))%>% select(gene, sample_id, disease_gene, exp.outlier) 
	require(plyr)
	test.out=rbind.fill(list(test.spli,test.exp))
	detach(package:plyr)
	test.out$outlier=pmin(rowSums(test.out[,c("spli.outlier", "exp.outlier")], na.rm=T),1)
	test.out=test.out %>% group_by(gene, disease_gene) %>% summarize(outlier=max(outlier)) %>% ungroup %>% select(disease_gene,outlier)
	

	test.glm=glm(test.out$disease_gene ~test.out$outlier, family=binomial)
	mod_summary =summary(test.glm)
	p.value=mod_summary[[12]][2, 4] # pvalue
	log_odds=mod_summary[[12]][2, 1] # coefficient estimate
	min.ci = mod_summary[[12]][2, 1] - (1.96*mod_summary[[12]][2, 2]) #
	max.ci = mod_summary[[12]][2, 1] + (1.96*mod_summary[[12]][2, 2]) #
	
	#summary.counts = as.data.frame(table(test.out))
    #if (nrow(summary.counts) != 4 | min(summary.counts$Freq) == 0) {
    #    cat("Warning: Skipping this threshold because zeros in contingency table.\n")
    #    return(NULL)
    #}
	#out.disease = summary.counts$Freq[summary.counts$disease_gene == 1 & summary.counts$outlier == 1]
	#nonout.disease = summary.counts$Freq[summary.counts$disease_gene == 1 & summary.counts$outlier == 0]
	#out.total = sum(summary.counts$Freq[summary.counts$outlier == 1])
	#nonout.total = sum(summary.counts$Freq[summary.counts$outlier == 0])
	#estimate = (out.disease/out.total)/(nonout.disease/nonout.total)
	## get bounds of confidence interval on the log of the proportion then exponentiate
	#log.se = sqrt(1/out.disease - 1/out.total + 1/nonout.disease - 1/nonout.total)
	#max.ci = estimate * exp(1.96*log.se)
	#min.ci = estimate * exp(-1.96*log.se)
	dfrow = list(LOG_ODDS=log_odds, CI.LOW=min.ci, CI.HIGH=max.ci, THRESH=thresh, COUNT = sum(test.out$outlier), PVALUE=p.value)
    return (dfrow)

}

proportion.ratios.glm.helper = function(counts, thresh) {
	test= counts %>% mutate(outlier=ifelse(abs(value)>=thresh, 1,0)) %>% select(gene, variable, disease_gene, outlier) %>% group_by(gene,variable) %>% summarize(disease_gene=max(disease_gene), outlier=max(outlier))%>% ungroup  %>% select(disease_gene,outlier)
	summary.counts = as.data.frame(table(test))
    if (nrow(summary.counts) != 4 | min(summary.counts$Freq) == 0) {
        cat("Warning: Skipping this threshold because zeros in contingency table.\n")
        return(NULL)
    }
	out.disease = summary.counts$Freq[summary.counts$disease_gene == 1 & summary.counts$outlier == 1]
	nonout.disease = summary.counts$Freq[summary.counts$disease_gene == 1 & summary.counts$outlier == 0]
	out.total = sum(summary.counts$Freq[summary.counts$outlier == 1])
	nonout.total = sum(summary.counts$Freq[summary.counts$outlier == 0])
	estimate = (out.disease/out.total)/(nonout.disease/nonout.total)
	# get bounds of confidence interval on the log of the proportion then exponentiate
	log.se = sqrt(1/out.disease - 1/out.total + 1/nonout.disease - 1/nonout.total)
	max.ci = estimate * exp(1.96*log.se)
	min.ci = estimate * exp(-1.96*log.se)
	dfrow = list(ESTIM=estimate, CI.LOW=min.ci, CI.HIGH=max.ci, THRESH=thresh, COUNT = sum(test$outlier))
    return (dfrow)
}


# generates HPO term list for each sample
get_list_hpo=function(sample, metadata){
	hpos=metadata %>% filter(sample_id==sample) %>% select(HPO_terms_ID)%>% pull
	hpos=unlist(strsplit(hpos,","))
	return(hpos)
}

# generates list of HPO alternative IDs from the input files
get_hpo_alt_id=function(myHPO_ID, altid_terms){
	if(altid_terms$ALT_ID[altid_terms$HPO_ID==myHPO_ID]!=""){
		alt_id=altid_terms %>% filter(HPO_ID==myHPO_ID) %>% mutate(ALT_ID=as.character(ALT_ID)) %>%pull
		hpos_alt=unlist(strsplit(alt_id, ","))
		hpos=c(myHPO_ID,hpos_alt)
	}else{hpos=myHPO_ID}
	return(hpos)

}


# function that look into list of alternative ids and expends original HPOs ids
get_alt_id=function(id,mylist){
if (id %in% names(mylist)){
	myhpos=c(id,mylist[[id]])

} else{
	myhpos=id
	for (element in names(mylist)){
		if (id %in% mylist[[element]]){
			myhpos=c(myhpos,mylist[[element]])
		}
	}
}
return(myhpos)
}

# function that look for all parent, child and alternative HPO terms from a given HPO term
process_sample_hpo=function(id){
	if(is.na(id) ){
		hpos=NA
		}
	else{ 
		parents=get_alt_id(id, parent_hpo_ids_list)
		childs=get_alt_id(id, child_hpo_ids_list)
		alt=get_alt_id(id,alt_hpo_ids_list)
		hpos=Reduce(union, list(parents,childs,alt))
		hpos=hpos[startsWith(as.character(hpos),"HP")]

}
	return(hpos)
}
# function that look for all HPO ters from a given samples and expend the list with parent, child and alt terms
process_sample_ids_hpo=function(sample, HPO_list){
	all_hpos=lapply(HPO_list[[sample]],process_sample_hpo)
	all_hpos=Reduce(union, all_hpos)
	return(all_hpos)
}
# function that generate a list of outlier genes with RV within 20bp of juction for which at list one HPO terms match with symptomes
get_hpo_comp=function(sample,outlier_RV, metadata, gene_to_pheno,pheno_annot, cadd_thres, HPOs_list){
	candidates=outlier_RV %>% filter(variable==sample) %>% filter(absZ>=2) %>% 
		#filter(distance_variant_junction<=20)%>%
		filter(as.numeric(as.character(phred_cadd))>=cadd_thres)%>%left_join(grch37, by=c("gene_name"="ensgene")) %>%
		select(gene_name, symbol, entrez, absZ, phred_cadd,gnomAD_AF_filt_withsgl)
hpo_comb=candidates %>% left_join(gene_to_pheno, by=c("entrez"="entrez")) %>% distinct %>% left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>% 
	mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPOs_list[[sample]], 1, 0)) %>% select(gene_name, symbol, absZ, phred_cadd,gnomAD_AF_filt_withsgl, hpo_in_input) %>% distinct %>%group_by(symbol, absZ, phred_cadd,gnomAD_AF_filt_withsgl)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>% filter(overlap_hpo>=1) %>% distinct%>% ungroup %>%distinct %>% group_by(symbol,overlap_hpo) %>% summarize(max_z=max(absZ)) %>% distinct %>% mutate(sample=sample)
return(as.data.frame(hpo_comb))
}
# function that generate a list of outlier genes with RV within 20bp of juction for which at list one HPO terms match with symptomes wondow approach
get_hpo_comp_window=function(sample,outlier_RV, metadata, gene_to_pheno,pheno_annot, cadd_thres, HPOs_list){
	if(cadd_thres==0){
	candidates=outlier_RV %>% filter(variable==sample) %>% filter(absZ>=2) %>% 
		filter(as.numeric(as.character(phred_cadd))>=cadd_thres | phred_cadd==".")%>%left_join(grch37, by=c("gene_name"="ensgene")) %>%
		select(gene_name, symbol, entrez, absZ, phred_cadd,gnomAD_AF_filt_withsgl)
hpo_comb=candidates %>% left_join(gene_to_pheno, by=c("entrez", "symbol")) %>% distinct %>% left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>% 
	mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPOs_list[[sample]], 1, 0)) %>% select(gene_name, symbol, absZ, phred_cadd,gnomAD_AF_filt_withsgl, hpo_in_input) %>% distinct %>%group_by(symbol, absZ, phred_cadd,gnomAD_AF_filt_withsgl)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>% filter(overlap_hpo>=1) %>% distinct%>% ungroup %>%distinct %>% group_by(symbol,overlap_hpo) %>% summarize(max_z=max(absZ)) %>% distinct %>% mutate(sample=sample)
}
else{
		candidates=outlier_RV %>% filter(variable==sample) %>% filter(absZ>=2) %>% 
		filter(as.numeric(as.character(phred_cadd))>=cadd_thres )%>%left_join(grch37, by=c("gene_name"="ensgene")) %>%
		select(gene_name, symbol, entrez, absZ, phred_cadd,gnomAD_AF_filt_withsgl)
hpo_comb=candidates %>% left_join(gene_to_pheno, by=c("entrez", "symbol")) %>% distinct %>% left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>% 
	mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPOs_list[[sample]], 1, 0)) %>% select(gene_name, symbol, absZ, phred_cadd,gnomAD_AF_filt_withsgl, hpo_in_input) %>% distinct %>%group_by(symbol, absZ, phred_cadd,gnomAD_AF_filt_withsgl)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>% filter(overlap_hpo>=1) %>% distinct%>% ungroup %>%distinct %>% group_by(symbol,overlap_hpo) %>% summarize(max_z=max(absZ)) %>% distinct %>% mutate(sample=sample)

}


return(as.data.frame(hpo_comb))
}


# function that generate a list of outlier genes for which at list one HPO terms match with symptomes
get_hpo_comp_nogen=function(sample,outlier, metadata, gene_to_pheno,pheno_annot, HPOs_list){
	candidates=outlier %>% filter(variable==sample) %>% filter(absZ>=2) %>% 
		left_join(grch37, by=c("gene"="ensgene")) %>%
		select(gene, symbol, entrez, absZ)
hpo_comb=candidates %>% left_join(gene_to_pheno, by=c("entrez", "symbol")) %>% distinct %>% left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>% 
	mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPOs_list[[sample]], 1, 0)) %>% select(gene, symbol, absZ, hpo_in_input) %>% distinct %>%group_by(symbol, absZ)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>% filter(overlap_hpo>=1) %>% distinct %>% ungroup %>%distinct %>% group_by(symbol,overlap_hpo) %>% summarize(max_z=max(absZ)) %>% distinct %>% mutate(sample=sample)
return(as.data.frame(hpo_comb))
}



# function that generate the number of outlier genes with RV within 20bp of juction for which at list one HPO terms match with symptomes. pLI score adjustable
get_hpo_match_nb<-function(sample,outlier_file,HPOs_list, pLI_thre){
 nb_match=outlier_file %>% filter(variable==sample)  %>% filter(absZ>=2) %>% 
 	left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct %>% 
 	filter(pLI>=pLI_thre) %>%
 	left_join(grch37, by=c("gene"="ensgene")) %>% select(gene, symbol, entrez, absZ) %>% 
 	left_join(gene_to_pheno, by=c("entrez", "symbol"))%>% distinct %>% 
 	left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>%
 	mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPOs_list[[sample]], 1, 0)) %>% 
 	select(gene, symbol, absZ, hpo_in_input) %>% distinct %>%group_by(gene,symbol, absZ)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>% 
 	filter(overlap_hpo>=1) %>% ungroup %>% unique%>% filter (! duplicated(gene)) %>% nrow 
return(nb_match)
}

# function that generate the number of outlier genes with RV within 20bp of juction for which at list one HPO terms match with symptomes, CADD and PLI scores adjustable
get_hpo_match_nb_RV=function(sample,outlier_RV, HPOs_list, cadd_thres, pLI_thre){
	nb_match=outlier_RV %>% filter(variable==sample) %>% filter(absZ>=2) %>% 
		filter(distance_variant_junction<=20)%>%
		filter(as.numeric(as.character(phred_cadd))>=cadd_thres)%>%
		left_join(exac,by=c("gene_name"="ensgene"), all.x=T) %>% distinct %>% 
		filter(pLI>=pLI_thre) %>%
		left_join(grch37, by=c("gene_name"="ensgene")) %>%
		select(gene, symbol, entrez, absZ) %>% 
		left_join(gene_to_pheno, by=c("entrez", "symbol"))%>% distinct %>% 
		left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>%
		mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPOs_list[[sample]], 1, 0)) %>%
		select(gene, symbol, absZ, hpo_in_input) %>% distinct %>%group_by(gene,symbol, absZ)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>% 
 		filter(overlap_hpo>=1) %>% ungroup %>% unique%>% filter (! duplicated(gene)) %>% nrow 
 	return(nb_match)
}

# function that generate the number of outlier genes with RV within 20bp of juction for which at list one HPO terms match with symptomes, CADD and PLI scores adjustable
get_hpo_match_nb_RV_window=function(sample,outlier_RV, HPOs_list, cadd_thres, pLI_thre){
	nb_match=outlier_RV %>% filter(variable==sample) %>% filter(absZ>=2) %>% 
		filter(as.numeric(as.character(phred_cadd))>=cadd_thres)%>%
		left_join(exac,by=c("gene_name"="ensgene"), all.x=T) %>% distinct %>% 
		filter(pLI>=pLI_thre) %>%
		left_join(grch37, by=c("gene_name"="ensgene")) %>%
		select(gene, symbol, entrez, absZ) %>% 
		left_join(gene_to_pheno, by=c("entrez", "symbol"))%>% distinct %>% 
		left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>%
		mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPOs_list[[sample]], 1, 0)) %>%
		select(gene, symbol, absZ, hpo_in_input) %>% distinct %>%group_by(gene,symbol, absZ)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>% 
 		filter(overlap_hpo>=1) %>% ungroup %>% unique%>% filter (! duplicated(gene)) %>% nrow 
 	return(nb_match)
}

data_summary <- function(x) {
   m <- mean(x)
   ymin <- m-sd(x)
   ymax <- m+sd(x)
   return(c(y=m,ymin=ymin,ymax=ymax))
}



# function that get the affected status of a sample
get_affected_status<-function(mysample, affected_status_df){
	affected_status=affected_status_df %>% filter(sample==mysample) %>% select(status)%>% pull
	return(affected_status)
}

# function that gets the institution where the smaple is from
get_institution<-function(mysample, metadata){
	institution=metadata %>% filter(sample_id==mysample) %>% select(institution)%>% pull
	return(institution)
}
# function that gets the kind of genetic data we have for the sample

get_technology<-function(mysample, metadata){
	technology=metadata %>% filter(sample_id==mysample) %>% select(variant_data)%>% pull
	return(technology)
}



# function that generate a list of outlier genes for which at list one HPO terms match with symptomes
get_hpo_comp_nogen=function(sample,outlier, metadata, gene_to_pheno,pheno_annot, HPOs_list){
	candidates=outlier %>% filter(variable==sample) %>% filter(absZ>=2) %>% 
		left_join(grch37, by=c("gene"="ensgene")) %>%
		select(gene, symbol, entrez, absZ)
hpo_comb=candidates %>% left_join(gene_to_pheno, by=c("entrez", "symbol")) %>% distinct %>% left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>% 
	mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPOs_list[[sample]], 1, 0)) %>% select(gene, symbol, absZ, hpo_in_input) %>% distinct %>%group_by(symbol, absZ)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>% filter(overlap_hpo>=1) %>% distinct %>% ungroup %>%distinct %>% group_by(symbol,overlap_hpo) %>% summarize(max_z=max(absZ)) %>% distinct %>% mutate(sample=sample)
return(as.data.frame(hpo_comb))
}

# function that generate a list of expression outlier genes for which at list one HPO terms match with symptomes
#exp_out_hpo=exp_outlier %>% filter(sample_id=="RD062")  %>% filter(abs(zscore)>=2) %>% left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct %>% filter(pLI>=0) %>%	left_join(grch37, by=c("gene"="ensgene")) %>% select(gene, symbol, entrez, zscore) %>% 	left_join(gene_to_pheno, by=c("entrez","symbol"))%>% distinct %>% 	left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>%	mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPO_extended[["RD062"]], 1, 0)) %>%select(gene, symbol, zscore, hpo_in_input) %>% distinct %>%group_by(gene,symbol, zscore)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>% filter(overlap_hpo>=1) %>% ungroup %>% unique%>% filter (! duplicated(gene)) %>% nrow,

get_hpo_comp_nogen_exp=function(sample,outlier, metadata, gene_to_pheno,pheno_annot, HPOs_list){
	candidates=outlier %>% filter(sample_id==sample) %>% filter(abs(zscore)>=2) %>% 
		left_join(grch37, by=c("gene"="ensgene")) %>%
		select(gene, symbol, entrez, zscore)
hpo_comb=candidates %>% left_join(gene_to_pheno, by=c("entrez", "symbol")) %>% distinct %>% left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>% 
	mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPOs_list[[sample]], 1, 0)) %>% select(gene, symbol, zscore, hpo_in_input) %>% distinct %>%group_by(symbol, zscore)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>% filter(overlap_hpo>=1) %>% distinct %>% ungroup %>%distinct %>% group_by(symbol,overlap_hpo) %>% summarize(max_z=max(abs(zscore))) %>% distinct %>% mutate(sample=sample)
return(as.data.frame(hpo_comb))
}


get_outlier_genes_exp<-function(mysample, threshold, outlier_df.m, pLI_thre, direction){
	if (direction=="over"){
	number_outliergenes=outlier_df.m %>% filter(sample_id==mysample)%>% filter(zscore>=threshold) %>%left_join(exac,by=c("gene"="ensgene"), all.x=T)%>% distinct %>% filter(pLI>=pLI_thre) %>%select(gene) %>% summarise(n_Elements = n_distinct(gene)) %>% pull
	}
	else if (direction=="under"){
	number_outliergenes=outlier_df.m %>% filter(sample_id==mysample)%>% filter(zscore<=-threshold) %>%left_join(exac,by=c("gene"="ensgene"), all.x=T)%>% distinct %>% filter(pLI>=pLI_thre) %>%select(gene) %>% summarise(n_Elements = n_distinct(gene)) %>% pull

	}
	return(number_outliergenes)

}

# function that generates the number of outlier genes for which we finnd ASE
get_outlier_ase_genes=function(mysample, outlier_df, threshold,gene_ase_list,pLI_thre, direction){
	if (direction=="over"){
	number_outliergenes=outlier_df %>% filter(sample_id==mysample)%>% filter(zscore>=threshold) %>%left_join(exac,by=c("gene"="ensgene"), all.x=T)%>% distinct %>% filter(pLI>=pLI_thre, gene %in% gene_ase_list[[mysample]]) %>%select(gene) %>% summarise(n_Elements = n_distinct(gene)) %>% pull
	}
	else if (direction=="under"){
	number_outliergenes=outlier_df%>% filter(sample_id==mysample)%>% filter(zscore<=-threshold) %>%left_join(exac,by=c("gene"="ensgene"), all.x=T)%>% distinct %>% filter(pLI>=pLI_thre, gene %in% gene_ase_list[[mysample]]) %>%select(gene) %>% summarise(n_Elements = n_distinct(gene)) %>% pull

	}
	return(number_outliergenes)

}

# function that generates the number of outlier genes for which we finnd ASE
get_outlier_river_genes=function(mysample, outlier_df, threshold,gene_river_list,pLI_thre, direction){
	if (direction=="over"){
	number_outliergenes=outlier_df %>% filter(sample_id==mysample)%>% filter(zscore>=threshold) %>%left_join(exac,by=c("gene"="ensgene"), all.x=T)%>% distinct %>% filter(pLI>=pLI_thre, gene %in% gene_river_list[[mysample]]) %>%select(gene) %>% summarise(n_Elements = n_distinct(gene)) %>% pull
	}
	else if (direction=="under"){
	number_outliergenes=outlier_df%>% filter(sample_id==mysample)%>% filter(zscore<=-threshold) %>%left_join(exac,by=c("gene"="ensgene"), all.x=T)%>% distinct %>% filter(pLI>=pLI_thre, gene %in% gene_river_list[[mysample]]) %>%select(gene) %>% summarise(n_Elements = n_distinct(gene)) %>% pull

	}
	return(number_outliergenes)

}

# function that generate the number of outlier genes with RV within 20bp of juction for which at list one HPO terms match with symptomes. pLI score adjustable
get_hpo_match_nb_exp<-function(sample,outlier_file,HPOs_list, pLI_thre, direction){
	if (direction=="over"){
 	nb_match=outlier_file %>% filter(sample_id==sample)  %>% filter(zscore>=2) %>% 
		left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct %>% 
		filter(pLI>=pLI_thre) %>%
		left_join(grch37, by=c("gene"="ensgene")) %>% select(gene, symbol, entrez, zscore) %>% 
		left_join(gene_to_pheno, by=c("entrez", "symbol"))%>% distinct %>% 
		left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>%
		mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPOs_list[[sample]], 1, 0)) %>% 
		select(gene, symbol, zscore, hpo_in_input) %>% distinct %>%group_by(gene,symbol, zscore)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>% 
		filter(overlap_hpo>=1) %>% ungroup %>% unique%>% filter (! duplicated(gene)) %>% nrow 
	}
	else if (direction=="under"){
	nb_match=outlier_file %>% filter(sample_id==sample)  %>% filter(zscore<=-2) %>% 
		left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct %>% 
		filter(pLI>=pLI_thre) %>%
		left_join(grch37, by=c("gene"="ensgene")) %>% select(gene, symbol, entrez, zscore) %>% 
		left_join(gene_to_pheno, by=c("entrez", "symbol"))%>% distinct %>% 
		left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>%
		mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPOs_list[[sample]], 1, 0)) %>% 
		select(gene, symbol, zscore, hpo_in_input) %>% distinct %>%group_by(gene,symbol, zscore)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>% 
		filter(overlap_hpo>=1) %>% ungroup %>% unique%>% filter (! duplicated(gene)) %>% nrow 
	}
	return(nb_match)
}

# function that generate the number of outlier genes with RV within 20bp of juction for which at list one HPO terms match with symptomes, CADD and PLI scores adjustable
get_hpo_match_nb_RV_exp=function(sample,outlier_RV, HPOs_list, cadd_thres, pLI_thre, direction){
	if (direction=="over"){
	nb_match=outlier_RV %>% filter(sample_id==sample) %>% filter(expressionZScore>=2) %>% 
		filter(variantDistanceFromGene<=10000)%>%
		filter(as.numeric(as.character(phred_cadd))>=cadd_thres)%>%
		left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct %>% 
		filter(pLI>=pLI_thre) %>%
		left_join(grch37, by=c("gene"="ensgene")) %>%
		select(gene, symbol, entrez, expressionZScore) %>% 
		left_join(gene_to_pheno, by=c("entrez", "symbol"))%>% distinct %>% 
		left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>%
		mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPOs_list[[sample]], 1, 0)) %>%
		select(gene, symbol, expressionZScore, hpo_in_input) %>% distinct %>%group_by(gene,symbol, expressionZScore)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>% 
 		filter(overlap_hpo>=1) %>% ungroup %>% unique%>% filter (! duplicated(gene)) %>% nrow 
	}
	if (direction=="under"){
	nb_match=outlier_RV %>% filter(sample_id==sample) %>% filter(expressionZScore<=-2) %>% 
		filter(variantDistanceFromGene<=10000)%>%
		filter(as.numeric(as.character(phred_cadd))>=cadd_thres)%>%
		left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct %>% 
		filter(pLI>=pLI_thre) %>%
		left_join(grch37, by=c("gene"="ensgene")) %>%
		select(gene, symbol, entrez, expressionZScore) %>% 
		left_join(gene_to_pheno, by=c("entrez", "symbol"))%>% distinct %>% 
		left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>%
		mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPOs_list[[sample]], 1, 0)) %>%
		select(gene, symbol, expressionZScore, hpo_in_input) %>% distinct %>%group_by(gene,symbol, expressionZScore)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>% 
 		filter(overlap_hpo>=1) %>% ungroup %>% unique%>% filter (! duplicated(gene)) %>% nrow 
	}
 	return(nb_match)
}

get_outlier_genes_RV_exp<-function(mysample,  outlier_df.m, distance,CADD, pLI_thre, direction){
	if (direction=="over"){
	number_outliergenes=outlier_df.m %>% filter(sample_id==mysample)%>% filter(expressionZScore>=2) %>% filter(variantDistanceFromGene<=distance)%>%
		 left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=CADD, pLI>=pLI_thre) %>% 
		select(gene) %>% summarise(n_Elements = n_distinct(gene)) %>% pull
	}
	if (direction=="under"){
		number_outliergenes=outlier_df.m %>% filter(sample_id==mysample)%>% filter(expressionZScore<=-2) %>% filter(variantDistanceFromGene<=distance)%>%
		 left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=CADD, pLI>=pLI_thre) %>% 
		select(gene) %>% summarise(n_Elements = n_distinct(gene)) %>% pull
	}
	return(number_outliergenes)

}

get_outlier_genes_RV_exp_10kb<-function(mysample,  outlier_df.m, CADD, pLI_thre, direction,myvariant_gene_pos){
	if (direction=="over"){
		if(CADD==0){
	number_outliergenes=outlier_df.m %>% filter(sample_id==mysample)%>% filter(expressionZScore>=2) %>% filter(variant_gene_pos==myvariant_gene_pos)%>%
		 left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=CADD | phred_cadd==".", pLI>=pLI_thre) %>% 
		select(gene) %>% summarise(n_Elements = n_distinct(gene)) %>% pull
	}else{
			number_outliergenes=outlier_df.m %>% filter(sample_id==mysample)%>% filter(expressionZScore>=2) %>% filter(variant_gene_pos==myvariant_gene_pos)%>%
		 left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=CADD, pLI>=pLI_thre) %>% 
		select(gene) %>% summarise(n_Elements = n_distinct(gene)) %>% pull

	}}

	if (direction=="under"){
		if(CADD==0){
		number_outliergenes=outlier_df.m %>% filter(sample_id==mysample)%>% filter(expressionZScore<=-2) %>% filter(variant_gene_pos==myvariant_gene_pos)%>%
		 left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=CADD|phred_cadd==".", pLI>=pLI_thre) %>% 
		select(gene) %>% summarise(n_Elements = n_distinct(gene)) %>% pull
	}else{
			number_outliergenes=outlier_df.m %>% filter(sample_id==mysample)%>% filter(expressionZScore<=-2) %>% filter(variant_gene_pos==myvariant_gene_pos)%>%
		 left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=CADD|phred_cadd==".", pLI>=pLI_thre) %>% 
		select(gene) %>% summarise(n_Elements = n_distinct(gene)) %>% pull
	
	}
	return(number_outliergenes)

}}

get_outlier_genes_RV_exp_10kb_ase<-function(mysample,  outlier_df.m, CADD, gene_ase_list,pLI_thre, direction,myvariant_gene_pos){
	if (direction=="over"){
			if(CADD==0){
	number_outliergenes=outlier_df.m %>% filter(sample_id==mysample)%>% filter(expressionZScore>=2) %>% filter(variant_gene_pos==myvariant_gene_pos)%>%
		 left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=CADD|phred_cadd==".", pLI>=pLI_thre, gene %in% gene_ase_list[[mysample]]) %>% 
		select(gene) %>% summarise(n_Elements = n_distinct(gene)) %>% pull
	}else{
	number_outliergenes=outlier_df.m %>% filter(sample_id==mysample)%>% filter(expressionZScore>=2) %>% filter(variant_gene_pos==myvariant_gene_pos)%>%
		 left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=CADD, pLI>=pLI_thre, gene %in% gene_ase_list[[mysample]]) %>% 
		select(gene) %>% summarise(n_Elements = n_distinct(gene)) %>% pull

		}}

	if (direction=="under"){
			if(CADD==0){
		number_outliergenes=outlier_df.m %>% filter(sample_id==mysample)%>% filter(expressionZScore<=-2) %>% filter(variant_gene_pos==myvariant_gene_pos)%>%
		 left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=CADD|phred_cadd==".", pLI>=pLI_thre, gene %in% gene_ase_list[[mysample]]) %>% 
		select(gene) %>% summarise(n_Elements = n_distinct(gene)) %>% pull
		}else{
	number_outliergenes=outlier_df.m %>% filter(sample_id==mysample)%>% filter(expressionZScore<=-2) %>% filter(variant_gene_pos==myvariant_gene_pos)%>%
	left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=CADD, pLI>=pLI_thre, gene %in% gene_ase_list[[mysample]]) %>% 
	select(gene) %>% summarise(n_Elements = n_distinct(gene)) %>% pull

		}
	}
	return(number_outliergenes)

}
get_outlier_genes_RV_exp_10kb_river<-function(mysample,  outlier_df.m, CADD, gene_river_list,pLI_thre, direction,myvariant_gene_pos){
	if (direction=="over"){
		if (CADD==0){
		number_outliergenes=outlier_df.m %>% filter(sample_id==mysample)%>% filter(expressionZScore>=2) %>% filter(variant_gene_pos==myvariant_gene_pos)%>%
		 left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=CADD|phred_cadd==".", pLI>=pLI_thre, gene %in% gene_river_list[[mysample]]) %>% 
		select(gene) %>% summarise(n_Elements = n_distinct(gene)) %>% pull
		
		}else{
		number_outliergenes=outlier_df.m %>% filter(sample_id==mysample)%>% filter(expressionZScore>=2) %>% filter(variant_gene_pos==myvariant_gene_pos)%>%
		 left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=CADD, pLI>=pLI_thre, gene %in% gene_river_list[[mysample]]) %>% 
		select(gene) %>% summarise(n_Elements = n_distinct(gene)) %>% pull

		}
	}
	if (direction=="under"){
		if(CADD==0){
					number_outliergenes=outlier_df.m %>% filter(sample_id==mysample)%>% filter(expressionZScore<=-2) %>% filter(variant_gene_pos==myvariant_gene_pos)%>%
		 left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=CADD|phred_cadd==".", pLI>=pLI_thre, gene %in% gene_river_list[[mysample]]) %>% 
		select(gene) %>% summarise(n_Elements = n_distinct(gene)) %>% pull

		}else{
			number_outliergenes=outlier_df.m %>% filter(sample_id==mysample)%>% filter(expressionZScore<=-2) %>% filter(variant_gene_pos==myvariant_gene_pos)%>%
			left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=CADD, pLI>=pLI_thre, gene %in% gene_river_list[[mysample]]) %>% 
			select(gene) %>% summarise(n_Elements = n_distinct(gene)) %>% pull

		}
	}
	return(number_outliergenes)

}

# function that generate the number of outlier genes with RV within 20bp of juction for which at list one HPO terms match with symptomes, CADD and PLI scores adjustable
get_hpo_match_nb_RV_exp_10kb=function(sample,outlier_RV, HPOs_list, cadd_thres, pLI_thre, direction,myvariant_gene_pos){
	if (direction=="over"){
		if(cadd_thres==0){
	nb_match=outlier_RV %>% filter(sample_id==sample) %>% filter(expressionZScore>=2) %>% 
		filter(variant_gene_pos==myvariant_gene_pos)%>%
		filter(as.numeric(as.character(phred_cadd))>=cadd_thres|phred_cadd==".")%>%
		left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct %>% 
		filter(pLI>=pLI_thre) %>%
		left_join(grch37, by=c("gene"="ensgene")) %>%
		select(gene, symbol, entrez, expressionZScore) %>% 
		left_join(gene_to_pheno, by=c("entrez", "symbol"))%>% distinct %>% 
		left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>%
		mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPOs_list[[sample]], 1, 0)) %>%
		select(gene, symbol, expressionZScore, hpo_in_input) %>% distinct %>%group_by(gene,symbol, expressionZScore)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>% 
 		filter(overlap_hpo>=1) %>% ungroup %>% unique%>% filter (! duplicated(gene)) %>% nrow 
	}else{
	nb_match=outlier_RV %>% filter(sample_id==sample) %>% filter(expressionZScore>=2) %>% 
		filter(variant_gene_pos==myvariant_gene_pos)%>%
		filter(as.numeric(as.character(phred_cadd))>=cadd_thres)%>%
		left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct %>% 
		filter(pLI>=pLI_thre) %>%
		left_join(grch37, by=c("gene"="ensgene")) %>%
		select(gene, symbol, entrez, expressionZScore) %>% 
		left_join(gene_to_pheno, by=c("entrez","symbol"))%>% distinct %>% 
		left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>%
		mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPOs_list[[sample]], 1, 0)) %>%
		select(gene, symbol, expressionZScore, hpo_in_input) %>% distinct %>%group_by(gene,symbol, expressionZScore)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>% 
 		filter(overlap_hpo>=1) %>% ungroup %>% unique%>% filter (! duplicated(gene)) %>% nrow 


	}
	}
	if (direction=="under"){
			if(cadd_thres==0){
	nb_match=outlier_RV %>% filter(sample_id==sample) %>% filter(expressionZScore<=-2) %>% 
		filter(variant_gene_pos==myvariant_gene_pos)%>%
		filter(as.numeric(as.character(phred_cadd))>=cadd_thres|phred_cadd==".")%>%
		left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct %>% 
		filter(pLI>=pLI_thre) %>%
		left_join(grch37, by=c("gene"="ensgene")) %>%
		select(gene, symbol, entrez, expressionZScore) %>% 
		left_join(gene_to_pheno, by=c("entrez","symbol"))%>% distinct %>% 
		left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>%
		mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPOs_list[[sample]], 1, 0)) %>%
		select(gene, symbol, expressionZScore, hpo_in_input) %>% distinct %>%group_by(gene,symbol, expressionZScore)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>% 
 		filter(overlap_hpo>=1) %>% ungroup %>% unique%>% filter (! duplicated(gene)) %>% nrow 
	}else{
	nb_match=outlier_RV %>% filter(sample_id==sample) %>% filter(expressionZScore<=-2) %>% 
		filter(variant_gene_pos==myvariant_gene_pos)%>%
		filter(as.numeric(as.character(phred_cadd))>=cadd_thres)%>%
		left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct %>% 
		filter(pLI>=pLI_thre) %>%
		left_join(grch37, by=c("gene"="ensgene")) %>%
		select(gene, symbol, entrez, expressionZScore) %>% 
		left_join(gene_to_pheno, by=c("entrez","symbol"))%>% distinct %>% 
		left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>%
		mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPOs_list[[sample]], 1, 0)) %>%
		select(gene, symbol, expressionZScore, hpo_in_input) %>% distinct %>%group_by(gene,symbol, expressionZScore)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>% 
 		filter(overlap_hpo>=1) %>% ungroup %>% unique%>% filter (! duplicated(gene)) %>% nrow 

	}}
 	return(nb_match)
}


# function that generate the number of outlier genes with RV within 20bp of juction for which at list one HPO terms match with symptomes, CADD and PLI scores adjustable
get_hpo_match_nb_RV_exp_10kb_ase=function(sample,outlier_RV, HPOs_list, gene_ase_list,cadd_thres, pLI_thre, direction,myvariant_gene_pos){
	if (direction=="over"){
		if(cadd_thres==0){
	nb_match=outlier_RV %>% filter(sample_id==sample) %>% filter(expressionZScore>=2) %>% 
		filter(variant_gene_pos==myvariant_gene_pos)%>%
		filter(as.numeric(as.character(phred_cadd))>=cadd_thres|phred_cadd==".")%>%
		left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct %>% 
		filter(pLI>=pLI_thre, gene %in% gene_ase_list[[sample]]) %>%
		left_join(grch37, by=c("gene"="ensgene")) %>%
		select(gene, symbol, entrez, expressionZScore) %>% 
		left_join(gene_to_pheno, by=c("entrez", "symbol"))%>% distinct %>% 
		left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>%
		mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPOs_list[[sample]], 1, 0)) %>%
		select(gene, symbol, expressionZScore, hpo_in_input) %>% distinct %>%group_by(gene,symbol, expressionZScore)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>% 
 		filter(overlap_hpo>=1) %>% ungroup %>% unique%>% filter (! duplicated(gene)) %>% nrow 
	}else{
	nb_match=outlier_RV %>% filter(sample_id==sample) %>% filter(expressionZScore>=2) %>% 
		filter(variant_gene_pos==myvariant_gene_pos)%>%
		filter(as.numeric(as.character(phred_cadd))>=cadd_thres)%>%
		left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct %>% 
		filter(pLI>=pLI_thre, gene %in% gene_ase_list[[sample]]) %>%
		left_join(grch37, by=c("gene"="ensgene")) %>%
		select(gene, symbol, entrez, expressionZScore) %>% 
		left_join(gene_to_pheno, by=c("entrez", "symbol"))%>% distinct %>% 
		left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>%
		mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPOs_list[[sample]], 1, 0)) %>%
		select(gene, symbol, expressionZScore, hpo_in_input) %>% distinct %>%group_by(gene,symbol, expressionZScore)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>% 
 		filter(overlap_hpo>=1) %>% ungroup %>% unique%>% filter (! duplicated(gene)) %>% nrow 

	}}
	if (direction=="under"){
		if(cadd_thres==0){
	nb_match=outlier_RV %>% filter(sample_id==sample) %>% filter(expressionZScore<=-2) %>% 
		filter(variant_gene_pos==myvariant_gene_pos)%>%
		filter(as.numeric(as.character(phred_cadd))>=cadd_thres|phred_cadd==".")%>%
		left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct %>% 
		filter(pLI>=pLI_thre, gene %in% gene_ase_list[[sample]]) %>%
		left_join(grch37, by=c("gene"="ensgene")) %>%
		select(gene, symbol, entrez, expressionZScore) %>% 
		left_join(gene_to_pheno, by=c("entrez", "symbol"))%>% distinct %>% 
		left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>%
		mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPOs_list[[sample]], 1, 0)) %>%
		select(gene, symbol, expressionZScore, hpo_in_input) %>% distinct %>%group_by(gene,symbol, expressionZScore)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>% 
 		filter(overlap_hpo>=1) %>% ungroup %>% unique%>% filter (! duplicated(gene)) %>% nrow 
	}else{
	nb_match=outlier_RV %>% filter(sample_id==sample) %>% filter(expressionZScore<=-2) %>% 
		filter(variant_gene_pos==myvariant_gene_pos)%>%
		filter(as.numeric(as.character(phred_cadd))>=cadd_thres)%>%
		left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct %>% 
		filter(pLI>=pLI_thre, gene %in% gene_ase_list[[sample]]) %>%
		left_join(grch37, by=c("gene"="ensgene")) %>%
		select(gene, symbol, entrez, expressionZScore) %>% 
		left_join(gene_to_pheno,by=c("entrez", "symbol"))%>% distinct %>% 
		left_join(pheno_annot,by=c("DiseaseId"="database_INDEX")) %>% distinct %>%
		mutate(hpo_in_input=ifelse(HPO_terms_ID %in% HPOs_list[[sample]], 1, 0)) %>%
		select(gene, symbol, expressionZScore, hpo_in_input) %>% distinct %>%group_by(gene,symbol, expressionZScore)%>% summarise(overlap_hpo=sum(hpo_in_input)) %>% 
 		filter(overlap_hpo>=1) %>% ungroup %>% unique%>% filter (! duplicated(gene)) %>% nrow 

	}}
 	return(nb_match)
}


process_variant=function(variant_freq){
	variant_freq=as.character(variant_freq)
	if (grepl(",", variant_freq)) {
		variant_freq=max(as.numeric(unlist(strsplit(variant_freq,","))))
		return(as.numeric(variant_freq))
	}
	else{ return(as.numeric(variant_freq))}
}

get_genes_with_highriver=function(sample,river_df){
	genes=river_df%>%filter(indID==sample)%>% filter(RIVERrank>=0.85)%>% select(Ensembl_id_DotStripped)%>% pull
	return(genes)
}

get_genes_with_ase=function(sample,ase_df){
	genes=ase_df%>%filter(sample_id==sample)%>% select(ensgene)%>% pull
	return(genes)
}
# generates list of HPO alternative IDs from the input files
get_hpo_alt_id=function(myHPO_ID, altid_terms){
	if(altid_terms$ALT_ID[altid_terms$HPO_ID==myHPO_ID]!=""){
		alt_id=altid_terms %>% filter(HPO_ID==myHPO_ID) %>% mutate(ALT_ID=as.character(ALT_ID)) %>%pull
		hpos_alt=unlist(strsplit(alt_id, ","))
		hpos=c(myHPO_ID,hpos_alt)
	}else{hpos=myHPO_ID}
	return(hpos)

}
# generates HPO term list for each sample
get_list_hpo=function(sample, metadata){
	hpos=metadata %>% filter(sample_id==sample) %>% select(HPO_terms_ID)%>% pull
	hpos=unlist(strsplit(hpos,","))
	return(hpos)
}
# function that look for all HPO ters from a given samples and expend the list with parent, child and alt terms
process_sample_ids_hpo=function(sample, HPO_list){
	all_hpos=lapply(HPO_list[[sample]],process_sample_hpo)
	all_hpos=Reduce(union, all_hpos)
	return(all_hpos)
}

# function that look for all parent, child and alternative HPO terms from a given HPO term
process_sample_hpo=function(id){
	if(is.na(id) ){
		hpos=NA
		}
	else{ 
		parents=get_alt_id(id, parent_hpo_ids_list)
		childs=get_alt_id(id, child_hpo_ids_list)
		alt=get_alt_id(id,alt_hpo_ids_list)
		hpos=Reduce(union, list(parents,childs,alt))
		hpos=hpos[startsWith(as.character(hpos),"HP")]

}
	return(hpos)
}
# function that look into list of alternative ids and expends original HPOs ids
get_alt_id=function(id,mylist){
if (id %in% names(mylist)){
	myhpos=c(id,mylist[[id]])

} else{
	myhpos=id
	for (element in names(mylist)){
		if (id %in% mylist[[element]]){
			myhpos=c(myhpos,mylist[[element]])
		}
	}
}
return(myhpos)
}



proportion.ratios.exp = function(zscore_table, thresholds) {
    # make empty data frame and fill it from rbinds (not very many, so it's fine)
    results = data.frame(ESTIM=numeric(), CI.LOW=numeric(), CI.HIGH=numeric(), THRESH = numeric(), COUNT = numeric(),
        stringsAsFactors = F)
    for (thresh in thresholds) {
        props = proportion.ratios.exp.helper (zscore_table, thresh)
        if (!is.null(props)) {
            results = rbind(results, props)
        }
    }
    return(results)
}


read_in_disease_lists<-function( disease_file,path_to_file){
	temp=read.table(paste0(path_to_file, disease_file), header=F)[,1]
	return(temp)
}


check_gene_names<-function(gene_list){
	if(all(startsWith(as.character(gene_list), "ENSG")) == TRUE){
		return(gene_list)
	}else {

		gene_list2=as.data.frame(gene_list) 
		names(gene_list2)=c('symbol')
		gene_list2=gene_list2 %>% merge(grch37,  by='symbol') %>% select(ensgene) %>% distinct
		gene_list2=gene_list2[,1]
		return(gene_list2)
	}
}



proportion.ratios.exp.helper = function(counts, thresh) {
	test= counts %>% mutate(outlier=ifelse(abs(zscore)>=thresh, 1,0)) %>% select(gene, disease_gene, outlier) 
	test.out=test %>% group_by(gene, disease_gene) %>% summarize(outlier=max(outlier)) %>% ungroup %>% select(disease_gene,outlier)
	test.glm=glm(test.out$disease_gene ~test.out$outlier, family=binomial)
	mod_summary =summary(test.glm)
	p.value=mod_summary[[12]][2, 4] # pvalue
	log_odds=mod_summary[[12]][2, 1] # coefficient estimate
	min.ci = mod_summary[[12]][2, 1] - (1.96*mod_summary[[12]][2, 2]) #
	max.ci = mod_summary[[12]][2, 1] + (1.96*mod_summary[[12]][2, 2]) #
	


	#summary.counts = as.data.frame(table(test))
    #if (nrow(summary.counts) != 4 | min(summary.counts$Freq) == 0) {
    #    cat("Warning: Skipping this threshold because zeros in contingency table.\n")
    #    return(NULL)
    #}
	#out.disease = summary.counts$Freq[summary.counts$disease_gene == 1 & summary.counts$outlier == 1]
	#nonout.disease = summary.counts$Freq[summary.counts$disease_gene == 1 & summary.counts$outlier == 0]
	#out.total = sum(summary.counts$Freq[summary.counts$outlier == 1])
	#nonout.total = sum(summary.counts$Freq[summary.counts$outlier == 0])
	#estimate = (out.disease/out.total)/(nonout.disease/nonout.total)
	## get bounds of confidence interval on the log of the proportion then exponentiate
	#log.se = sqrt(1/out.disease - 1/out.total + 1/nonout.disease - 1/nonout.total)
	#max.ci = estimate * exp(1.96*log.se)
	#min.ci = estimate * exp(-1.96*log.se)
	dfrow = list(LOG_ODDS=log_odds, CI.LOW=min.ci, CI.HIGH=max.ci, THRESH=thresh, COUNT = sum(test.out$outlier), PVALUE=p.value)
    return (dfrow)
}

emp_ast <- function(pvalues) {
	sapply(pvalues, function(i) {
		if (i > 0.1) return("")
		else if (i<=0.1 & i >0.05)return(".")
		else if (i <= 0.05 & i > 0.01) return("*")
		else if (i <= 0.01 & i > 0.001) return("**")
		else return("***")
		})
}


### for enrichement HPO analysis

get_outlier_genes_names<-function(mysample, threshold, outlier_df.m, pLI_thre){
	outliergenes=outlier_df.m %>% filter(variable==mysample)%>% filter(absZ>=threshold) %>%left_join(exac,by=c("gene"="ensgene"), all.x=T)%>% distinct %>% filter(pLI>=pLI_thre) %>%select(gene) %>% filter (! duplicated(gene)) %>% pull
	return(outliergenes)
}

get_outlier_genes_RV_window_names<-function(mysample, threshold, outlier_df.m,CADD, pLI_thre){
	if(CADD==0){
		outliergenes=outlier_df.m %>% filter(variable==mysample)%>% filter(absZ>=threshold) %>%
		left_join(exac,by=c("gene_name"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=CADD| phred_cadd==".", pLI>=pLI_thre) %>% 
		select(gene_name) %>%  filter (! duplicated(gene_name)) %>% pull
	}
	else{
		outliergenes=outlier_df.m %>% filter(variable==mysample)%>% filter(absZ>=threshold) %>%
		left_join(exac,by=c("gene_name"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=CADD, pLI>=pLI_thre) %>% 
		select(gene_name) %>% filter (! duplicated(gene_name)) %>% pull

	}
	return(outliergenes)
}



is_hpo=function(sample_id, gene){
	if (sample_id %in% HPO_extended.df$sample){
		temp=HPO_extended.df %>% filter(sample==sample_id)
		hpo=ifelse(gene %in% temp$ensgene, 1, 0)
		return(hpo)

}else{return(NA)}
}

get_outlier_genes_RV_exp_names<-function(mysample,  outlier_df.m, distance,CADD, pLI_thre, direction){
	if (direction=="over"){
	outliergenes=outlier_df.m %>% filter(sample_id==mysample)%>% filter(expressionZScore>=2) %>% filter(variantDistanceFromGene<=distance)%>%
		 left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=CADD, pLI>=pLI_thre) %>% 
		select(gene)  %>% filter (! duplicated(gene)) %>% pull
	}
	if (direction=="under"){
		outliergenes=outlier_df.m %>% filter(sample_id==mysample)%>% filter(expressionZScore<=-2) %>% filter(variantDistanceFromGene<=distance)%>%
		 left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=CADD, pLI>=pLI_thre) %>% 
		select(gene) %>% filter (! duplicated(gene)) %>% pull
	}
	return(outliergenes)

}

get_outlier_genes_RV_exp_10kb_names<-function(mysample,  outlier_df.m, CADD, pLI_thre, direction,myvariant_gene_pos){
	if (direction=="over"){
		if(CADD==0){
	outliergenes=outlier_df.m %>% filter(sample_id==mysample)%>% filter(expressionZScore>=2) %>% filter(variant_gene_pos==myvariant_gene_pos)%>%
		 left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=CADD | phred_cadd==".", pLI>=pLI_thre) %>% 
		select(gene) %>% filter (! duplicated(gene)) %>% pull
	}else{
		outliergenes=outlier_df.m %>% filter(sample_id==mysample)%>% filter(expressionZScore>=2) %>% filter(variant_gene_pos==myvariant_gene_pos)%>%
		left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=CADD, pLI>=pLI_thre) %>% 
		select(gene) %>% filter (! duplicated(gene)) %>% pull

	}}

	if (direction=="under"){
		if(CADD==0){
		outliergenes=outlier_df.m %>% filter(sample_id==mysample)%>% filter(expressionZScore<=-2) %>% filter(variant_gene_pos==myvariant_gene_pos)%>%
		 left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=CADD|phred_cadd==".", pLI>=pLI_thre) %>% 
		select(gene) %>% filter (! duplicated(gene)) %>% pull
	}else{
			outliergenes=outlier_df.m %>% filter(sample_id==mysample)%>% filter(expressionZScore<=-2) %>% filter(variant_gene_pos==myvariant_gene_pos)%>%
		 left_join(exac,by=c("gene"="ensgene"), all.x=T) %>% distinct%>%filter(as.numeric(as.character(phred_cadd))>=CADD|phred_cadd==".", pLI>=pLI_thre) %>% 
		select(gene) %>% filter (! duplicated(gene)) %>% pull
	
	}

}
	return(outliergenes)
}



get_outlier_genes_exp_names<-function(mysample, threshold, outlier_df.m, pLI_thre, direction){
	if (direction=="over"){
	outliergenes=outlier_df.m %>% filter(sample_id==mysample)%>% filter(zscore>=threshold) %>%left_join(exac,by=c("gene"="ensgene"), all.x=T)%>% distinct %>% filter(pLI>=pLI_thre) %>%select(gene) %>%  filter (! duplicated(gene)) %>% pull
	}
	else if (direction=="under"){
	outliergenes=outlier_df.m %>% filter(sample_id==mysample)%>% filter(zscore<=-threshold) %>%left_join(exac,by=c("gene"="ensgene"), all.x=T)%>% distinct %>% filter(pLI>=pLI_thre) %>%select(gene) %>%  filter (! duplicated(gene)) %>% pull

	}
	return(outliergenes)

}

get_percent_HPO_genes_splicing_RV=function(samples_set, filter_function, Zscore_threshold, outlier_file,CADD_threshold,pLI_threshold, method_description){
	gene_list=lapply(samples_set, filter_function, threshold=Zscore_threshold, outlier_df.m=outlier_file, CADD=CADD_threshold, pLI_thre=pLI_threshold)
	names(gene_list)=samples_set
	gene_list=compact(gene_list)
	gene_list <- Map(cbind, gene_list, sample = names(gene_list))
	gene_list.df=as.data.frame(do.call("rbind", gene_list))
	gene_list.df$hpo=mapply(is_hpo,gene_list.df$sample,gene_list.df$V1)
	gene_list_hpo_percent=gene_list.df %>% group_by(sample) %>% mutate(n_genes=n()) %>% mutate(percent_hpo=sum(hpo)/n_genes)%>% select(sample,n_genes,percent_hpo)%>% filter(!duplicated(sample)) 
	gene_list_hpo_percent$method=method_description
return(gene_list_hpo_percent)
}


get_percent_HPO_genes_expression=function(samples_set, filter_function, Zscore_threshold, outlier_file,pLI_threshold, out_direction,method_description){
	gene_list=lapply(samples_set, filter_function, threshold=Zscore_threshold, outlier_df.m=outlier_file,  pLI_thre=pLI_threshold, direction=out_direction)
	names(gene_list)=samples_set
	gene_list=compact(gene_list)
	gene_list <- Map(cbind, gene_list, sample = names(gene_list))
	gene_list.df=as.data.frame(do.call("rbind", gene_list))
	gene_list.df$hpo=mapply(is_hpo,gene_list.df$sample,gene_list.df$V1)
	gene_list_hpo_percent=gene_list.df %>% group_by(sample) %>% mutate(n_genes=n()) %>% mutate(percent_hpo=sum(hpo)/n_genes)%>% select(sample,n_genes,percent_hpo)%>% filter(!duplicated(sample)) 
	gene_list_hpo_percent$method=method_description
	gene_list_hpo_percent$direction=out_direction
	gene_list_hpo_percent$variant_position=NA

return(gene_list_hpo_percent)
}

get_percent_HPO_genes_expression_RV=function(samples_set, filter_function, outlier_file,CADD_threshold,pLI_threshold, out_direction,method_description, variant_pos){
	gene_list=lapply(samples_set, filter_function,  outlier_df.m=outlier_file,  CADD=CADD_threshold, pLI_thre=pLI_threshold, direction=out_direction, myvariant_gene_pos=variant_pos)
	names(gene_list)=samples_set
	gene_list=compact(gene_list)
	gene_list <- Map(cbind, gene_list, sample = names(gene_list))
	gene_list.df=as.data.frame(do.call("rbind", gene_list))
	gene_list.df$hpo=mapply(is_hpo,gene_list.df$sample,gene_list.df$V1)
	gene_list_hpo_percent=gene_list.df %>% group_by(sample) %>% mutate(n_genes=n()) %>% mutate(percent_hpo=sum(hpo)/n_genes)%>% select(sample,n_genes,percent_hpo)%>% filter(!duplicated(sample)) 
	gene_list_hpo_percent$method=method_description
	gene_list_hpo_percent$direction=out_direction
	gene_list_hpo_percent$variant_position=variant_pos


return(gene_list_hpo_percent)
}

get_percent=function(values, threshold){
	return(length(which(values<threshold))/length(values)*100)
}

emp_ast <- function(pvalues) {
	sapply(pvalues, function(i) {
		if (i > 0.05) return("")
		else if (i <= 0.5 & i > 0.01) return("*")
		else if (i <= 0.01 & i > 0.001) return("**")
		else return("***")
		})
}

get_plot_labels <- function(n) {
	iterator <- c(0.5, 0.25, 0.1, 10/n, 5/n, 2/n, 1/n)
	nsamples <- round(n*iterator)
	nsamples <- c(nsamples[order(nsamples)], nsamples[-1])
	iterator <- unique(c(iterator, 1-iterator))
	iterator <- iterator[order(iterator)]
	for_return <- matrix(nrow=2, ncol=length(iterator))
	for_return[1, ] <- iterator
	for_return[2, ] <- nsamples
	return(for_return)
}

ggplot_build <- function(dat, nsamples, iterator, emp) {
	p1 <- ggplot(dat) + 
	geom_hline(yintercept=0) + 
	geom_errorbar(aes(x=threshold, ymin=error_low, ymax=error_high), colour="gray", size=0.5, width=0.5) + 
	geom_point(aes(x=threshold, y=real_coefficient, fill=predictor), shape=21, colour="Black", size=2.5) +
	facet_grid(predictor ~ ., scales="free_y") +
	labs(x="Percentile", y="Log Odds Ratio") +
	scale_x_continuous(breaks=seq(1, length(nsamples), 1), labels=paste0(round(iterator*100, 3), "% (n=", nsamples, ")")) +
	geom_text(aes(x=threshold, y=real_coefficient+0.01), label=emp, size=3) +
	theme_bw() +
	theme(strip.background=element_blank(),
		panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank(),
		axis.text=element_text(size=9),
		axis.text.x=element_text(angle=45, hjust=1),
		panel.border=element_blank()) +
	annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
	annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) +
	guides(fill=FALSE) +
	scale_fill_manual(values=c("indianred3", "royalblue3", "orange2")) 
	return(p1)
}

ggplot_build_z <- function(dat, iterator, emp, ngenes) {
	p1 <- ggplot(dat) + 
	geom_hline(yintercept=0) + 
	geom_errorbar(aes(x=threshold, ymin=error_low, ymax=error_high), colour="gray", size=0.5, width=0.5) + 
	geom_point(aes(x=threshold, y=real_coefficient, fill=predictor), shape=21, colour="Black", size=2.5) +
	facet_grid(predictor ~ ., scales="free_y") +
	labs(x="Percentile", y="Log Odds Ratio") +
	scale_x_continuous(breaks=seq(1, length(iterator), 1), labels=iterator) +
	geom_text(aes(x=threshold, y=real_coefficient+0.1), label=emp, size=3) +
	geom_text(aes(x=threshold, y=error_high+0.2), label=ngenes, size=3, angle=45) +
	theme_bw() +
	theme(strip.background=element_blank(),
		panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank(),
		axis.text=element_text(size=9),
		panel.border=element_blank()) +
	annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1) +
	annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1) +
	guides(fill=FALSE) +
	scale_fill_manual(values=c("indianred3", "royalblue3", "orange2")) 
	return(p1)
}

# Process junction files for figure s2
process_junc_file=function(sample, junction_dir,junc_suffix){
	print(sample)
	temp_junc=paste0(junction_dir,sample,junc_suffix)
	temp_junc=read_tsv(temp_junc, col_names=F) 
	colnames(temp_junc)=c("chr", "intron_start", "intron_end", "strand", "intron_motif", "annotation_status", "uniquely_mapped", "multi_mapping", "max_alignment_overhang", "gene_chr", "gene_start", "gene_end", "gene_name")
	temp_junc=temp_junc %>% filter(annotation_status==1)%>% mutate(junction=paste(chr, intron_start,intron_end,gene_name, sep="_"))%>% select(junction, uniquely_mapped) 
	return(temp_junc)
}
