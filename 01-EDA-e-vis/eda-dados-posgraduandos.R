# Você precisará instalar esses pacotes. Faça install.packages("nome") para cada um.
library(dplyr, warn.conflicts = F)
library(readr)
library(ggplot2)
# theme_set(theme_bw()) # você pode preferir os gráficos assim
library(gmodels)

# ====================================
# LER, ARRUMAR, LIMPAR
# ====================================

# Repare que estou usando readr::read_csv em vez de base::read.csv
# read_csv é bem mais rápido, adivinha melhor os tipos das colunas
# e nunca usa Factor, sempre usa String no lugar. Geralmente ajuda.
dados = read_csv("dados//Dados de alunos para as aulas de FPCC-report.csv")

View(dados)
# usando dplyr
dados %>% View()

# Renomeia as colunas e mantém apenas as que quero
dados = dados %>% 
  select(curso = `De que curso você é aluno?`, 
         sexo = `Você é...`, 
         altura = `Qual a sua altura em centímetros?`,
         repositorios = `Em quantos repositórios de software você lembra ter contribuído nos últimos 2 anos?`, 
         linguagens = `Em quantas linguagens de programação você se considera fluente?`,
         projetos_de_pesquisa = `Em quantos projetos de pesquisa você lembra ter participado?`, 
         confianca_estatistica = `Qual seu nível de confiança hoje no uso de métodos estatísticos para analisar o resultado de um experimento?`, 
         gosta_de_forms = `O quanto você gosta de formulários online? (Obrigado por ter respondido este!)`,
         submissao = `Submit Date (UTC)`,
         fpcc2 = `Você já cursou, está cursando ou não cursou FPCC 2?`)

# Remove NAs
dados = dados %>% 
  filter(complete.cases(dados))

str(dados)

# ====================================
# EXPLORAR, VISUALIZAR
# ====================================

# ------------------------------------
# Variáveis numéricas
# ------------------------------------

# Altura
ggplot(dados, mapping = aes(x = "altura", y = altura)) + 
  geom_point()

ggplot(dados, aes(x = "altura", y = altura)) + 
  geom_violin() + 
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 3, alpha = 0.5)

ggplot(dados, mapping = aes(x = altura)) + 
  #geom_histogram(bins = 12) 
  # geom_freqpoly(bins = 20)
  geom_density() 
  # geom_rug()

# Linguagens de programação
dados %>% 
  ggplot(mapping = aes(x = "Quantas", y = linguagens)) + 
  # geom_point()
  geom_count()

dados %>% 
  ggplot(mapping = aes(x = linguagens)) + 
  geom_histogram(bins = 6)

# Compare com repos
dados %>% 
  #filter(repositorios < 10) %>% 
  ggplot(mapping = aes(x = repositorios)) + 
  geom_histogram(bins = 5) + 
  geom_vline(xintercept = mean(dados$repositorios), colour = "orange")

# Qual o formato esperado da distribuição para as demais variáveis?

# Validade e confiabilidade das variáveis

# Medias e medianas

dados %>% 
  #filter(repositorios < 10) %>% 
  ggplot(mapping = aes(x = repositorios, fill = curso)) + 
  geom_histogram(bins = 6) + 
  facet_grid(curso ~ .) + 
  geom_rug()

dados %>% 
  group_by(curso) %>% 
  summarise(linguagens.media = mean(linguagens), 
            repos.medio = mean(repositorios))

dados %>% 
  group_by(sexo) %>% 
  summarise(altura.media = mean(altura), 
            sd.altura = sd(altura))


# ------------------------------------
# CATEGÓRICAS
# ------------------------------------

# Curso e sexo são categóricos
ggplot(dados) + 
  geom_bar(mapping = aes(x = curso), stat = "count") + 
  coord_flip()

ggplot(dados) + 
  geom_bar(mapping = aes(x = sexo), stat = "count") + 
  coord_flip()

ggplot(dados) + 
  geom_bar(mapping = aes(x = curso, fill = sexo), position = "fill") + # tente position = "stack"/"dodge"/"fill"
  coord_flip()


# ------------------------------------
# DUAS VARIÁVEIS
# ------------------------------------
ggplot(dados, aes (y = repositorios, 
                   x = linguagens)) + 
  geom_point(size = 2, alpha = 0.6)

ggplot(dados, aes(x = curso, y = altura)) + 
  #geom_boxplot(alpha = 0.2) +
  geom_violin() + 
  geom_point(position = position_jitter(width = 0.07), size = 4, alpha = 0.5)

CrossTable(dados$sexo, dados$curso, prop.chisq = FALSE)

