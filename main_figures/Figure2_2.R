
fsize=15
RD_theme=theme_classic()+
	theme(axis.text.x= element_text(size=fsize),
	axis.text.y= element_text(size=fsize), 
	axis.title = element_text(size = fsize), 
    legend.text = element_text(size = fsize-1), 
    legend.title = element_text(size = fsize), 
    axis.ticks = element_line(size = 0.1),
	strip.background=element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	axis.text=element_text(size=9),
	panel.border=element_blank()) 


# plot results
filter_withsglt.exp.case.10kb.plot=ggplot(sample_exp_outlier_filter_withsgl_10kb_up.df.m.cases, aes(x=variable, y=value, fill=variable) )+ 
	geom_boxplot(color="black", notch=F, show.legend = FALSE)+
	geom_point(size = 0, stroke = 0)+
	scale_fill_manual(values=rep("lightgrey", 13), breaks=c("prop_pLI","prop_RIVER","prop_ASE", "prop_HPO","prop_RV", "prop_RV.pLI", "prop_CADD","prop_RV.RIVER","prop_RV.HPO", "prop_RV.ASE", "prop_RV.HPO.CADD","prop_RV.HPO.CADD.pLI","prop_RV.HPO.CADD.ASE"),
		labels=c(expression("pLI " >="0.9"), expression("RIVER " >= "0.85"), "ASE", "HPO match", expression("Rare variant within 10kb","1 + 5","5 + CADD score  " >= "10"), "2 + 5", "4 + 5", "3 + 5", "4 + 7", "1 + 2 + 5 + 7", "2 + 4 + 5 + 7"), name="Filter")+
	scale_x_discrete(breaks=c("prop_pLI","prop_RIVER","prop_ASE", "prop_HPO","prop_RV", "prop_RV.pLI", "prop_CADD","prop_RV.RIVER","prop_RV.HPO", "prop_RV.ASE", "prop_RV.HPO.CADD","prop_RV.HPO.CADD.pLI","prop_RV.HPO.CADD.ASE"),
		labels=c("1", "2", "3","4","5", "6", "7", "8", "9", "10", "11", "12", "13"))+
	labs(x="Filter", y="Proportion of outlier genes") + RD_theme +guides(fill = guide_legend(override.aes = list(size = 4,shape =c(49,50,51,52,53,54,55,56,57,58,59,60,61))))+theme(legend.position = c(0.8,0.7))
filter_withsglt.exp.case.10kb.plot

ggsave('fig_2_prop_candidates_filters_10kb_upstream.pdf', filter_withsglt.exp.case.10kb.plot, path='/srv/scratch/restricted/rare_diseases/analysis/manuscript/figures_revision/', width=6, height=6)


filter_withsglt.exp.case.10kb_techno.plot=ggplot(sample_exp_outlier_filter_withsgl_10kb_up.df.m.cases, aes(x=variable, y=value, fill=variable, color=technology) )+ 
	geom_boxplot( notch=F, show.legend = FALSE)+
	geom_point(size = 0, stroke = 0)+
	scale_fill_manual(values=rep("lightgrey", 13), breaks=c("prop_pLI","prop_RIVER","prop_ASE", "prop_HPO","prop_RV", "prop_RV.pLI", "prop_CADD","prop_RV.RIVER","prop_RV.HPO", "prop_RV.ASE", "prop_RV.HPO.CADD","prop_RV.HPO.CADD.pLI","prop_RV.HPO.CADD.ASE"),
		labels=c(expression("pLI " >="0.9"), expression("RIVER " >= "0.85"), "ASE", "HPO match", expression("Rare variant within 10kb","1 + 5","5 + CADD score  " >= "10"), "2 + 5", "4 + 5", "3 + 5", "4 + 7", "1 + 11", "3 + 12"), name="Filter")+
	scale_color_manual(values=c("#003366", "#CCCC99"))+ 
	scale_x_discrete(breaks=c("prop_pLI","prop_RIVER","prop_ASE", "prop_HPO","prop_RV", "prop_RV.pLI", "prop_CADD","prop_RV.RIVER","prop_RV.HPO", "prop_RV.ASE", "prop_RV.HPO.CADD","prop_RV.HPO.CADD.pLI","prop_RV.HPO.CADD.ASE"),
		labels=c("1", "2", "3","4","5", "6", "7", "8", "9", "10", "11", "12", "13"))+
	labs(x="Filter", y="Proportion of outlier genes") + RD_theme +guides(fill = guide_legend(override.aes = list(size = 4,shape =c(49,50,51,52,53,54,55,56,57,58,59,60,61))), colour=guide_legend(override.aes = list(size = 4,shape =c(1,1)), fill=c("#003366", "#CCCC99")))+theme(legend.position = c(0.8,0.5))
filter_withsglt.exp.case.10kb_techno.plot

ggsave('sup_prop_candidates_filters_10kb_upstream_techno.pdf', filter_withsglt.exp.case.10kb_techno.plot, path='/srv/scratch/restricted/rare_diseases/analysis/manuscript/figures_revision/', width=6, height=6)

filter_withsglt.exp.case.10kb_affected_status.plot=ggplot(sample_exp_outlier_filter_withsgl_10kb_up.df.m, aes(x=variable, y=value, fill=variable, color=affected_status) )+ 
	geom_boxplot( notch=F, show.legend = FALSE)+
	geom_point(size = 0, stroke = 0)+
	scale_fill_manual(values=rep("lightgrey", 13), breaks=c("prop_pLI","prop_RIVER","prop_ASE", "prop_HPO","prop_RV", "prop_RV.pLI", "prop_CADD","prop_RV.RIVER","prop_RV.HPO", "prop_RV.ASE", "prop_RV.HPO.CADD","prop_RV.HPO.CADD.pLI","prop_RV.HPO.CADD.ASE"),
		labels=c(expression("pLI " >="0.9"), expression("RIVER " >= "0.85"), "ASE", "HPO match", expression("Rare variant within 10kb","1 + 5","5 + CADD score  " >= "10"), "2 + 5", "4 + 5", "3 + 5", "4 + 7", "1 + 11", "3 + 12"), name="Filter")+
	scale_color_brewer(palette="Paired",name="Affected status")+ 
	scale_x_discrete(breaks=c("prop_pLI","prop_RIVER","prop_ASE", "prop_HPO","prop_RV", "prop_RV.pLI", "prop_CADD","prop_RV.RIVER","prop_RV.HPO", "prop_RV.ASE", "prop_RV.HPO.CADD","prop_RV.HPO.CADD.pLI","prop_RV.HPO.CADD.ASE"),
		labels=c("1", "2", "3","4","5", "6", "7", "8", "9", "10", "11", "12", "13"))+
	labs(x="Filter", y="Proportion of outlier genes") + RD_theme +guides(fill = guide_legend(override.aes = list(size = 4,shape =c(49,50,51,52,53,54,55,56,57,58,59,60,61))), colour=guide_legend(override.aes = list(size = 4,shape =c(1,1)), fill=c("#003366", "#CCCC99")))+theme(legend.position = c(0.8,0.5))
filter_withsglt.exp.case.10kb_affected_status.plot

ggsave('sup_prop_candidates_filters_10kb_upstream_affected_status.pdf', filter_withsglt.exp.case.10kb_affected_status.plot, path='/srv/scratch/restricted/rare_diseases/analysis/manuscript/figures_revision/', width=6, height=6)



color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
mycol=c("steelblue","dodgerblue1","goldenrod2","wheat2","moccasin","thistle1", 
"lightgoldenrod4","lavenderblush", "orange1",
"cornsilk2", "mediumblue","darkseagreen1",
"tomato","mediumpurple4", "lightskyblue2",
"khaki1","sienna2", "honeydew4",
"dodgerblue2","cyan3", "cyan4",
"lightgoldenrod3","darkseagreen4", "darkorange1",  
"olivedrab1","bisque1", "orangered1",
"green1","ivory3","purple3",
"skyblue", "deepskyblue1",  "springgreen",  
"orchid2", "oldlace", "bisque3",
"hotpink4","lightcoral","lightsteelblue2", 
"tomato3", "seagreen4", "wheat1", 
"khaki", "lightsalmon4",  "bisque4",
"darkseagreen2", "hotpink1","darkgoldenrod2",
"pink3", "coral2","palevioletred2",
"lightsalmon1",  "violetred3","cornsilk3",
"maroon3", "dodgerblue","tan3", 
"honeydew3", "seashell3", "seashell1",
"navy",  "darkorchid3","darkolivegreen",
"lightsteelblue","darksalmon","lightgoldenrod",
"orangered4","cornflowerblue","honeydew", 
"khaki4","peachpuff4","red2", 
"pink1", "steelblue1","rosybrown3",
"lightblue", "lightpink3","turquoise4",
"firebrick1","aquamarine","lavenderblush4",
"palevioletred1","antiquewhite",  "blue", 
"coral1","navajowhite2",  "mediumspringgreen",
"red3",  "indianred3","salmon", 
"lightsteelblue3","mintcream", "chartreuse","gold3")



all_ind_plot= ggplot(exp_outlier_number.df, aes(filter,value))+ 
 	geom_quasirandom(aes(y =value, x = filter, colour = sample))+
 	stat_summary(fun.data = median.quartile, geom = "pointrange",position=position_nudge(x=0.5,y=0))+
	coord_flip()+
	#geom_path(group=sample, colour="grey")+
 	theme(legend.position="")+
 	scale_colour_manual(values=mycol) +
 	labs(y="log10(Number of genes +1)",x="Filter")+
 	scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) 


ggsave('Fig_2D_expression_outlier_all_ind_plot.pdf', all_ind_plot, path='/srv/scratch/restricted/rare_diseases/analysis/manuscript/figures_revision/', width=7, height=7)
