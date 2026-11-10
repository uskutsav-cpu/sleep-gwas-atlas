load('work/sleep_pfi_charls/CHARLS_2011_2018_PFI_Sleep_Data.rda')
for (v in c('age','BMI','score','Glucose','Creatinine','Urea','TG','HDL','LDL')) g[[v]] <- as.numeric(as.character(g[[v]]))
f <- PFI ~ Qnight + age + sex + residence + BMI + education + smoke + drink + diabetes + dyslipidemia + hypertension + Glucose + Creatinine + Urea + TG + HDL + LDL
m <- glm(f, data=g, family=binomial())
s <- summary(m)$coefficients[c('QnightQ1','QnightQ3'),,drop=FALSE]
out <- data.frame(contrast=rownames(s), n=nobs(m), events=sum(m$y), odds_ratio=exp(s[,1]), lower95=exp(s[,1]-1.96*s[,2]), upper95=exp(s[,1]+1.96*s[,2]), p=s[,4], row.names=NULL)
write.csv(out,'outputs/stage6_charls_pfi_available_case_sensitivity.csv',row.names=FALSE)
print(out)
