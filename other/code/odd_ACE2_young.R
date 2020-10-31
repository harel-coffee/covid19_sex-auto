expr_meta_subset = expr_meta[!is.na(expr_meta$age) &  expr_meta$age == "0-19" & expr_meta$ACE2 > 5 ,]

expr_meta_subset = expr_meta[expr_meta$ESR1 ,]

source = names(table(expr_meta$source_name_ch1)[table(expr_meta$source_name_ch1) > 30])
expr_meta_subset = expr_meta[expr_meta$source_name_ch1 %in% source, ]
tail(sort(by(expr_meta_subset$ACE2, expr_meta_subset$source_name_ch1 , median)), 50)

tail(sort(by(expr_meta$ACE2, expr_meta$source_name , median)), 50)
