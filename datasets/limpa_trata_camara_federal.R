library(readr)
library(dplyr)
detalhes = read_csv("camara_federal/deputados-detalhes.csv")
votacoes = read_csv("camara_federal/votacoes.csv") %>%
    select(-cunha)
resumo_votacoes = votacoes %>% 
    select(1:9) %>% 
    distinct() %>% 
    write_csv("camara_federal/resumo_votacoes.csv")
votacoes %>% 
    select(-5, -8) %>% 
    write_csv("camara_federal/votacoes-menor.csv")
proposicoes = read_csv("camara_federal/proposicoes.csv")

gastos = rbind(read_csv("../../gastos-parlamentares/2015-2016-ate-mes-7.csv"), 
      read_csv("../../gastos-parlamentares/resultado.csv"), 
      read_csv("../../gastos-parlamentares/resultado-2017.csv"))
gastos %>% distinct() %>% write_csv("camara_federal/gastos-cota_atividade_parlamentar.csv")
