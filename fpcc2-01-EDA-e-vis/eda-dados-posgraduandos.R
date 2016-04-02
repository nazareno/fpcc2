require(ggplot2)
require(dplyr)

dados <- read.csv("dados//Dados de alunos para as aulas de FPCC-report.csv")

dados <- filter(dados, complete.cases(dados))
str(dados)

ggplot(dados, aes(x = "altura", y = Qual.a.sua.altura.em.centímetros.)) + 
  geom_violin() + 
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 3, alpha = 0.5)

ggplot(dados, aes(x = De.que.curso.você.é.aluno.)) + 
         geom_histogram()

ggplot(dados, aes(x = Você.é...)) + 
  geom_histogram()

ggplot(dados, aes (y = Em.quantos.repositórios.de.software.você.lembra.ter.contribuído.nos.últimos.2.anos., 
                   x = Em.quantas.linguagens.de.programação.você.se.considera.fluente.)) + 
  geom_point(size = 2, alpha = 0.6)



ggplot(dados, aes(x = De.que.curso.você.é.aluno., y = Qual.a.sua.altura.em.centímetros.)) + 
  #geom_boxplot(alpha = 0.2) +
  geom_violin() + 
  geom_point(position = position_jitter(width = 0.07), size = 4, alpha = 0.5)

require(gmodels)
CrossTable(dados$Você.é..., dados$De.que.curso.você.é.aluno., prop.chisq = FALSE)

ggplot(dados, aes(x = Qual.a.sua.altura.em.centímetros.)) + 
  geom_density() + 
  geom_rug(alpha = 0.8) + 
  geom_vline(xintercept = mean(dados$Qual.a.sua.altura.em.centímetros.), color = "Blue") + 
  geom_vline(xintercept = median(dados$Qual.a.sua.altura.em.centímetros.), color = "Darkorange") + 
  theme_bw()

ggplot(dados, aes(x = Em.quantos.repositórios.de.software.você.lembra.ter.contribuído.nos.últimos.2.anos.)) + 
  geom_bar(binwidth = 1) + 
  geom_vline(xintercept = mean(dados$Em.quantos.repositórios.de.software.você.lembra.ter.contribuído.nos.últimos.2.anos.), color = "Blue") + 
  geom_vline(xintercept = median(dados$Em.quantos.repositórios.de.software.você.lembra.ter.contribuído.nos.últimos.2.anos.), color = "Darkorange") + 
  theme_bw()
