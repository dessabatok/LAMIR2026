# Instalação, caso necessário
install.packages(c("FactoMineR", "factoextra", "dlm", "dplyr", "ggplot2", "Metrics"))

# Carregamento
library(FactoMineR)
library(factoextra)
library(dlm)
library(dplyr)
library(ggplot2)
library(Metrics)
library(tidyr)

# Carregar a base de dados
df = read.csv("C:/Users/karol/OneDrive/Área de Trabalho/Mestrado/Códigos/base_completa.csv", sep = ";")

# música 1

feat1 <- df %>% 
  filter(file==2) %>% 
  select(-frameTime, -Arousal, -Valence, -file)

pca1 <- prcomp(feat1, center = TRUE, scale. = TRUE, rank. = 10)
summary(pca1)
feat1pca <- pca1$x

a1 <- df %>% filter(file == 2) %>% pull(Arousal)
v1 <- df %>% filter(file == 2) %>% pull(Valence)
Y1  <- cbind(Arousal = a1, Valence = v1)   # matriz n x 2, limpa

X1 <- cbind(Intercept = 1, as.matrix(feat1pca))

n <- nrow(Y1)
p <- ncol(X1)
m <- 2  # dimensão da observação: Arousal e Valence

# FF bloco-diagonal: coeficientes próprios para cada série
FF  <- matrix(0, nrow = m, ncol = 2*p)
JFF <- matrix(0, nrow = m, ncol = 2*p)
JFF[1, 1:p]         <- 1:p   # equação do Arousal usa as colunas de X1
JFF[2, (p+1):(2*p)] <- 1:p   # equação do Valence usa as MESMAS colunas de X1 (reaproveitadas)
GG <- diag(2*p)

# --- função de construção do modelo ---
buildFun <- function(parm) {
  # V (2x2) via Cholesky -> garante positiva definida
  L <- matrix(c(exp(parm[1]), 0,
                parm[2],      exp(parm[3])), nrow = 2, byrow = TRUE)
  V <- L %*% t(L)
  
  # W diagonal (2p variâncias, log-parametrizadas p/ positividade)
  Wdiag <- exp(parm[4:(3 + 2*p)])
  
  dlm(
    FF = FF, JFF = JFF, X = X1,
    GG = GG,
    V  = V,
    W  = diag(Wdiag, nrow = 2*p),
    m0 = rep(0, 2*p),
    C0 = diag(1e7, 2*p)   # prior difusa
  )
}

# --- estimação ---
parm_init <- rep(0, 3 + 2*p)  # 3 p/ V (Cholesky) + 2p p/ W

fit <- dlmMLE(Y1, parm = parm_init, build = buildFun, hessian = TRUE)
fit$convergence  # 0 = convergiu; != 0 = revisar chute inicial / maxit

modelo_bivariado <- buildFun(fit$par)

# --- filtro de Kalman ---
filtro <- dlmFilter(Y1, modelo_bivariado)

# --- MSE (equivalente ao getMSE() do pydlm: erro de previsão 1 passo à frente) ---
burn_in <- 5  # ajuste conforme quantas observações o gráfico mostrar como instabilidade

mse_arousal <- mean((Y1[-(1:burn_in), 1] - filtro$f[-(1:burn_in), 1])^2, na.rm = TRUE)
mse_valence <- mean((Y1[-(1:burn_in), 2] - filtro$f[-(1:burn_in), 2])^2, na.rm = TRUE)

suave  <- dlmSmooth(filtro)

# --- 1) Observado vs. previsão 1 passo à frente ---
df_pred <- data.frame(
  tempo = 1:n,
  Arousal_obs  = Y1[,1],
  Arousal_prev = filtro$f[,1],
  Valence_obs  = Y1[,2],
  Valence_prev = filtro$f[,2]
) %>%
  pivot_longer(-tempo, names_to = c("serie", "tipo"), names_sep = "_", values_to = "valor")

ggplot(df_pred, aes(x = tempo, y = valor, color = tipo)) +
  geom_line() +
  facet_wrap(~serie, ncol = 1, scales = "free_y") +
  labs(title = "Observado vs. previsão 1 passo à frente",
       x = "Tempo", y = NULL, color = NULL) +
  theme_minimal()



####### Caso onde W é não correlacionada #######

# FF/JFF com estados intercalados: (coef1_A, coef1_V, coef2_A, coef2_V, ...)
FF  <- matrix(0, nrow = 2, ncol = 2*p)
JFF <- matrix(0, nrow = 2, ncol = 2*p)
for (i in 1:p) {
  JFF[1, 2*i - 1] <- i   # equação Arousal usa a coluna i de X1
  JFF[2, 2*i]     <- i   # equação Valence usa a MESMA coluna i de X1
}
GG <- diag(2*p)

buildFun <- function(parm) {
  # V (2x2): covariância residual entre Arousal e Valence
  L_V <- matrix(c(exp(parm[1]), 0,
                  parm[2],      exp(parm[3])), nrow = 2, byrow = TRUE)
  V <- L_V %*% t(L_V)
  
  # W: p blocos 2x2 -- cada bloco correlaciona a evolução do MESMO
  # coeficiente entre as duas equações (ex: intercepto_A com intercepto_V)
  idx <- 4
  W_blocks <- vector("list", p)
  for (i in 1:p) {
    a <- parm[idx]; b <- parm[idx + 1]; c <- parm[idx + 2]
    L_i <- matrix(c(exp(a), 0, b, exp(c)), nrow = 2, byrow = TRUE)
    W_blocks[[i]] <- L_i %*% t(L_i)
    idx <- idx + 3
  }
  W <- as.matrix(bdiag(W_blocks))
  
  dlm(FF = FF, JFF = JFF, X = X1, GG = GG, V = V, W = W,
      m0 = rep(0, 2*p), C0 = diag(1e7, 2*p))
}

parm_init <- rep(0, 3 + 3*p)  # 3 p/ V + 3 p/ bloco de W, vezes p blocos

fit <- dlmMLE(Y1, parm = parm_init, build = buildFun, hessian = TRUE,
              control = list(maxit = 500))
fit$convergence

modelo_bivariado <- buildFun(fit$par)
# --- filtro de Kalman ---
filtro <- dlmFilter(Y1, modelo_bivariado)

# --- MSE (equivalente ao getMSE() do pydlm: erro de previsão 1 passo à frente) ---
mse_arousal <- mean((Y1[,1] - filtro$f[,1])^2, na.rm = TRUE)
mse_valence <- mean((Y1[,2] - filtro$f[,2])^2, na.rm = TRUE)

filtro <- dlmFilter(Y1, modelo_bivariado)
suave  <- dlmSmooth(filtro)

# --- 1) Observado vs. previsão 1 passo à frente ---
df_pred <- data.frame(
  tempo = 1:n,
  Arousal_obs  = Y1[,1],
  Arousal_prev = filtro$f[,1],
  Valence_obs  = Y1[,2],
  Valence_prev = filtro$f[,2]
) %>%
  pivot_longer(-tempo, names_to = c("serie", "tipo"), names_sep = "_", values_to = "valor")

ggplot(df_pred, aes(x = tempo, y = valor, color = tipo)) +
  geom_line() +
  facet_wrap(~serie, ncol = 1, scales = "free_y") +
  labs(title = "Observado vs. previsão 1 passo à frente",
       x = "Tempo", y = NULL, color = NULL) +
  theme_minimal()

