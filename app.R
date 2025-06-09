# ============================================================================
# APLICAÇÃO SHINY PARA ANÁLISE DE EXPRESSÃO GÉNICA DIFERENCIAL COM DESeq2
# Análise de dados RNA-seq de células MCF7 tratadas com Mefloquina e Withaferina A
# ============================================================================

# Carregamento das bibliotecas necessárias
library(shiny)
library(shinydashboard)
library(DT)
library(DESeq2)
library(ggplot2)
library(pheatmap)  
library(RColorBrewer)
library(plotly)
library(dplyr)
library(data.table)
library(ggrepel)
library(factoextra)
#library(vsn)


# ============================================================================
# PROCESSAMENTO DOS DADOS EM BACKGROUND (SERVIDOR) - VERSÃO DEBUG
# ============================================================================

prepare_dataset <- function() {
  
  cat("=== INICIANDO PREPARAÇÃO DOS DADOS ===\n")
  
  # Verificar se os ficheiros existem
  required_files <- c("sample_info.txt", "gene_count_matrix.csv")
  
  for(file in required_files) {
    if(!file.exists(file)) {
      stop(paste("Ficheiro não encontrado:", file, 
                 "\nCertifique-se de que todos os ficheiros estão na pasta da aplicação."))
    }
  }
  
  # Carregar informação das amostras
  cat("1. Carregando informação das amostras...\n")
  coldata <- read.delim("sample_info.txt", sep = "\t", stringsAsFactors = FALSE)
  
  # Debug: mostrar estrutura inicial
  cat("   - Dimensões coldata inicial:", dim(coldata), "\n")
  cat("   - Nomes das colunas:", paste(colnames(coldata), collapse = ", "), "\n")
  cat("   - Primeiras 3 linhas:\n")
  print(head(coldata, 3))
  
  # Definir row names corretamente
  if(ncol(coldata) > 1) {
    row.names(coldata) <- coldata[, 1]
    coldata <- coldata[, -1, drop = FALSE]
    cat("   - Row names definidos, nova dimensão:", dim(coldata), "\n")
  } else {
    stop("Erro: sample_info.txt deve ter pelo menos 2 colunas")
  }
  
  # Filtrar para grupo 5 (MW - Mefloquina e Withaferina A)
  cat("2. Aplicando filtro de grupo...\n")
  group <- 5
  total_samples <- nrow(coldata)
  expected_samples <- 9 * 6  # 9 amostras por grupo, 6 grupos
  
  cat("   - Total de amostras disponíveis:", total_samples, "\n")
  cat("   - Amostras esperadas para 6 grupos:", expected_samples, "\n")
  
  if(total_samples < expected_samples) {
    cat("   - AVISO: Número de amostras menor que esperado. Usando todas as amostras disponíveis.\n")
    # Usar todas as amostras se não houver amostras suficientes
  } else {
    # Filtrar conforme especificado
    selected_indices <- c(1:9, (9*group+1):(9*group+9))
    cat("   - Índices selecionados:", paste(selected_indices, collapse = ", "), "\n")
    selected_indices <- selected_indices[selected_indices <= total_samples]
    coldata <- coldata[selected_indices, , drop = FALSE]
    cat("   - Amostras após filtro:", nrow(coldata), "\n")
  }
  
  # Carregar matriz de contagens
  cat("3. Carregando matriz de contagens...\n")
  if(file.exists("gene_count_matrix.csv")) {
    cts <- read.csv("gene_count_matrix.csv", row.names = 1, check.names = FALSE)
    cat("   - Usando gene_count_matrix.csv\n")
  } else if(file.exists("transcript_count_matrix.csv")) {
    cts <- read.csv("transcript_count_matrix.csv", row.names = 1, check.names = FALSE)
    cat("   - Usando transcript_count_matrix.csv\n")
    
    # Se usar transcript_count_matrix, processar anotações
    if(file.exists("transcripts_annot.txt")) {
      cat("   - Processando anotações de transcritos...\n")
      annot <- read.delim("transcripts_annot.txt")
      annot$new <- paste0(annot$transcript_id, "|", annot$gene_name)
      cts <- merge(annot, cts, by.x = "transcript_id", by.y = 0, all.y = FALSE)
      row.names(cts) <- cts[, 3]
      cts <- cts[, -c(1, 2, 3)]
    }
  } else {
    stop("Nenhum ficheiro de matriz de contagens encontrado!")
  }
  
  cat("   - Dimensões cts inicial:", dim(cts), "\n")
  cat("   - Primeiros 5 nomes de colunas:", paste(head(colnames(cts), 5), collapse = ", "), "\n")
  
  # Limpeza dos nomes das colunas
  cat("4. Limpando nomes das colunas...\n")
  colnames_original <- colnames(cts)
  colnames(cts) <- gsub("\\.gtf_", "", colnames(cts))
  colnames(cts) <- gsub("\\.bam$", "", colnames(cts))
  colnames(cts) <- gsub("^X", "", colnames(cts))
  
  if(!identical(colnames_original, colnames(cts))) {
    cat("   - Nomes alterados. Primeiros 5 novos nomes:", paste(head(colnames(cts), 5), collapse = ", "), "\n")
  } else {
    cat("   - Nenhuma alteração nos nomes das colunas\n")
  }
  
  # Verificar sobreposição entre nomes
  cat("5. Verificando correspondência de nomes...\n")
  cat("   - Amostras em coldata:", nrow(coldata), "\n")
  cat("   - Amostras em cts:", ncol(cts), "\n")
  cat("   - Primeiros 5 nomes em coldata:", paste(head(row.names(coldata), 5), collapse = ", "), "\n")
  cat("   - Primeiros 5 nomes em cts:", paste(head(colnames(cts), 5), collapse = ", "), "\n")
  
  common_samples <- intersect(colnames(cts), row.names(coldata))
  cat("   - Amostras em comum:", length(common_samples), "\n")
  
  if(length(common_samples) == 0) {
    cat("   - ERRO: Nenhuma amostra em comum encontrada!\n")
    cat("   - Tentando correspondência aproximada...\n")
    
    for(i in 1:min(5, length(row.names(coldata)))) {
      sample_name <- row.names(coldata)[i]
      matches <- grep(sample_name, colnames(cts), value = TRUE)
      if(length(matches) > 0) {
        cat("     - Amostra", sample_name, "pode corresponder a:", paste(matches, collapse = ", "), "\n")
      }
    }
    
    stop("Erro: Nenhuma correspondência encontrada entre amostras")
  } else {
    cat("   - Amostras em comum:", paste(head(common_samples, 5), collapse = ", "), "...\n")
  }
  
  # Filtrar para amostras comuns
  cat("6. Filtrando para amostras comuns...\n")
  cts <- cts[, colnames(cts) %in% row.names(coldata)]
  coldata <- coldata[row.names(coldata) %in% colnames(cts), , drop = FALSE]
  
  cat("   - Dimensões após filtro:\n")
  cat("     - cts:", dim(cts), "\n")
  cat("     - coldata:", dim(coldata), "\n")
  
  # Filtros de qualidade
  cat("7. Aplicando filtros de qualidade...\n")
  
  # Remover amostras com baixa cobertura (reduzido para 1M reads)
  low_coverage_cut <- 1000000
  sample_totals <- colSums(cts)
  cat("   - Totais de reads por amostra (primeiros 5):", paste(head(sample_totals, 5), collapse = ", "), "\n")
  
  low_coverage_samples <- names(sample_totals[sample_totals < low_coverage_cut])
  
  if(length(low_coverage_samples) > 0) {
    cat("   - Removendo", length(low_coverage_samples), "amostras com baixa cobertura\n")
    cat("   - Amostras removidas:", paste(low_coverage_samples, collapse = ", "), "\n")
    cts <- cts[, !colnames(cts) %in% low_coverage_samples]
    coldata <- coldata[!row.names(coldata) %in% low_coverage_samples, , drop = FALSE]
  } else {
    cat("   - Nenhuma amostra removida por baixa cobertura\n")
  }
  
  # Remover genes com contagens muito baixas
  initial_genes <- nrow(cts)
  cts <- cts[rowSums(cts) > 10, ]  # Pelo menos 10 reads no total
  cts <- cts[rowSums(cts > 0) >= 3, ]  # Presente em pelo menos 3 amostras
  
  cat("   - Genes removidos por baixa expressão:", initial_genes - nrow(cts), "\n")
  cat("   - Genes restantes:", nrow(cts), "\n")
  
  # Ordenar para garantir correspondência exata
  cat("8. Ordenando para correspondência final...\n")
  common_samples_final <- intersect(colnames(cts), row.names(coldata))
  cts <- cts[, common_samples_final]
  coldata <- coldata[common_samples_final, , drop = FALSE]
  
  # Verificação final
  cat("9. Verificação final:\n")
  cat("   - ncol(cts):", ncol(cts), "\n")
  cat("   - nrow(coldata):", nrow(coldata), "\n")
  cat("   - Nomes correspondem exatamente:", all(colnames(cts) == row.names(coldata)), "\n")
  
  if(!all(colnames(cts) == row.names(coldata))) {
    cat("   - ERRO: Correspondência final falhou!\n")
    cat("   - Primeiros 5 nomes cts:", paste(head(colnames(cts), 5), collapse = ", "), "\n")
    cat("   - Primeiros 5 nomes coldata:", paste(head(row.names(coldata), 5), collapse = ", "), "\n")
    stop("ERRO FINAL: Correspondência entre amostras falhou!")
  }
  
  # Preparar fatores para análise
  cat("10. Preparando fatores para análise...\n")
  variable_of_interest <- 3  # Treatment_TimePoint column
  
  if(ncol(coldata) >= variable_of_interest) {
    cat("   - Usando coluna", variable_of_interest, ":", colnames(coldata)[variable_of_interest], "\n")
    
    # Verificar valores únicos na coluna
    unique_values <- unique(coldata[, variable_of_interest])
    cat("   - Valores únicos:", paste(unique_values, collapse = ", "), "\n")
    
    # Definir fatores específicos para o estudo
    expected_levels <- c("DMSO_0", "DMSO_9", "DMSO_24", "MW_0", "MW_9", "MW_24")
    present_levels <- intersect(expected_levels, unique_values)
    
    if(length(present_levels) > 0) {
      coldata[, variable_of_interest] <- factor(coldata[, variable_of_interest], levels = present_levels)
      cat("   - Níveis do fator definidos:", paste(present_levels, collapse = ", "), "\n")
    } else {
      # Usar todos os valores únicos como fatores
      coldata[, variable_of_interest] <- as.factor(coldata[, variable_of_interest])
      cat("   - Usando todos os valores únicos como fatores\n")
    }
    
    cat("   - Distribuição das amostras:\n")
    print(table(coldata[, variable_of_interest]))
    
  } else {
    cat("   - AVISO: Coluna", variable_of_interest, "não encontrada. Usando primeira coluna.\n")
    variable_of_interest <- 1
    coldata[, variable_of_interest] <- as.factor(coldata[, variable_of_interest])
  }
  
  cat("=== PREPARAÇÃO DOS DADOS CONCLUÍDA COM SUCESSO! ===\n")
  cat("Dimensões finais: cts =", dim(cts), ", coldata =", dim(coldata), "\n")
  
  return(list(cts = cts, coldata = coldata, variable_of_interest = variable_of_interest))
}

# Função para executar análise DESeq2 com debug
run_deseq2_analysis <- function(cts, coldata, variable_of_interest, reference = NULL, experimental = NULL) {
  
  cat("=== INICIANDO ANÁLISE DESeq2 ===\n")
  
  # Verificação final antes de criar DESeqDataSet
  cat("1. Verificação final dos dados:\n")
  cat("   - ncol(cts):", ncol(cts), "\n")
  cat("   - nrow(coldata):", nrow(coldata), "\n")
  cat("   - Correspondência de nomes:", all(colnames(cts) == rownames(coldata)), "\n")
  
  if(ncol(cts) != nrow(coldata)) {
    stop(paste("Erro de dimensões: ncol(cts) =", ncol(cts), 
               "!= nrow(coldata) =", nrow(coldata)))
  }
  
  if(!all(colnames(cts) == rownames(coldata))) {
    stop("Erro: nomes das amostras não correspondem exatamente")
  }
  
  # Definir referência se fornecida
  if(!is.null(reference)) {
    cat("2. Definindo nível de referência:", reference, "\n")
    coldata[, variable_of_interest] <- relevel(coldata[, variable_of_interest], reference)
  }
  
  # Criar design formula
  design_formula <- as.formula(paste("~", colnames(coldata)[variable_of_interest]))
  cat("3. Fórmula do design:", deparse(design_formula), "\n")
  
  # Criar DESeqDataSet
  cat("4. Criando DESeqDataSet...\n")
  tryCatch({
    dds <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = coldata,
                                  design = design_formula)
    
    cat("   - DESeqDataSet criado com sucesso!\n")
    cat("   - Amostras:", ncol(dds), "\n")
    cat("   - Genes:", nrow(dds), "\n")
    
    # Executar análise
    cat("5. Executando análise DESeq2...\n")
    dds <- DESeq(dds)
    
    cat("6. Obtendo resultados...\n")
    # Obter nomes das comparações disponíveis
    comparison_names <- resultsNames(dds)
    cat("   - Comparações disponíveis:", paste(comparison_names, collapse = ", "), "\n")
    
    # Se experimental especificado, encontrar comparação correspondente
    if(!is.null(experimental)) {
      which_comparison <- grep(experimental, comparison_names)
      if(length(which_comparison) > 0) {
        cat("   - Usando comparação:", comparison_names[which_comparison[1]], "\n")
        res <- results(dds, 
                       name = comparison_names[which_comparison[1]],
                       pAdjustMethod = "fdr",
                       lfcThreshold = 1,
                       alpha = 0.05)
        comparison_name <- comparison_names[which_comparison[1]]
      } else {
        cat("   - Comparação experimental não encontrada, usando primeira disponível\n")
        res <- results(dds, alpha = 0.05)
        comparison_name <- comparison_names[1]
      }
    } else {
      # Usar primeira comparação disponível
      res <- results(dds, alpha = 0.05)
      comparison_name <- comparison_names[1]
    }
    
    # Transformações para visualização
    cat("7. Aplicando transformações para visualização...\n")
    vsd <- vst(dds, blind = FALSE)
    
    cat("=== ANÁLISE DESeq2 CONCLUÍDA COM SUCESSO! ===\n")
    cat("Genes diferencialmente expressos (p.adj < 0.05):", sum(res$padj < 0.05, na.rm = TRUE), "\n")
    
    return(list(dds = dds, res = res, vsd = vsd, comparison_name = comparison_name))
    
  }, error = function(e) {
    cat("ERRO na análise DESeq2:\n")
    cat("Mensagem:", e$message, "\n")
    stop(e)
  })
}

# ============================================================================
# INTERFACE DO UTILIZADOR (UI)
# ============================================================================

ui <- dashboardPage(
  skin = "blue",
  
  # Cabeçalho
  dashboardHeader(
    title = tags$div(
      style = "display: flex; align-items: center; height: 50px;",
      tags$a(href = "#", class = "sidebar-toggle", `data-toggle` = "offcanvas", role = "button",
             style = "display: flex; align-items: center; justify-content: center; width: 40px; height: 40px; margin-right: 10px;",
             icon("fa-bars")
      ),
      tags$span("RNA-seq: MCF7 MW", 
                style = "margin-left: 10px; font-weight: bold; font-size: 18px; line-height: 1; display: flex; align-items: center;")
    ),
    titleWidth = 300
  ),
  
  # Barra lateral
  dashboardSidebar(
    width = 300,
    sidebarMenu(
      menuItem("Informações", tabName = "methods", icon = icon("info-circle")),
      menuItem("Análise Principal", tabName = "main_analysis", icon = icon("chart-line")),
      menuItem("Heatmaps", tabName = "heatmaps", icon = icon("fire"))
    ),
    br(),
    h4("Parâmetros", style = "color: white; margin-left: 15px;"),
    
    selectInput("comparison_type", 
                "Escolher Comparação:",
                choices = list(
                  "MW_24h vs MW_0h (Efeito temporal)" = "MW_24_vs_MW_0",
                  "MW_24h vs DMSO_24h (Efeito tratamento)" = "MW_24_vs_DMSO_24",
                  "MW_9h vs DMSO_9h (Efeito tratamento 9h)" = "MW_9_vs_DMSO_9",
                  "MW_0h vs DMSO_0h (Controlo basal)" = "MW_0_vs_DMSO_0"
                ),
                selected = "MW_24_vs_MW_0"),
    
    
    
    
    
    numericInput("pval_threshold", 
                 "Threshold p-value ajustado:", 
                 value = 0.05, 
                 min = 0.001, 
                 max = 0.1, 
                 step = 0.001),
    
    numericInput("lfc_threshold", 
                 "Threshold Log2FoldChange:", 
                 value = 1, 
                 min = 0, 
                 max = 3, 
                 step = 0.1),
    
    
    numericInput("n_genes_heatmap", 
                 "Número de genes no heatmap:", 
                 value = 25, 
                 min = 10, 
                 max = 50, 
                 step = 5),
    
    # Botão de exportação
    downloadButton("download_results", 
                   "Exportar Resultados", 
                   class = "btn-primary",
                   style = "width: 90%; margin-left: 15px;")
  ),
  
  # Corpo principal
  dashboardBody(
    tags$head(
      tags$style(HTML("
        /* Esconder o toggle original da navbar */
        .main-header .navbar .sidebar-toggle {
          display: none !important;
        }
        
        /* Garantir que o header tenha altura fixa */
        .main-header {
          height: 50px !important;
        }
        
        /* Centrar verticalmente todo o conteúdo do logo */
        .main-header .logo {
          display: flex !important;
          align-items: center !important;
          height: 50px !important;
          padding: 0 15px !important;
        }
        
        /* Estilizar o ícone das barras */
        .sidebar-toggle {
          color: white !important;
          text-decoration: none !important;
          padding: 0 !important;
          border-radius: 4px;
          transition: background-color 0.3s ease;
        }
        
        .sidebar-toggle:hover {
          background-color: rgba(255, 255, 255, 0.1) !important;
          color: white !important;
        }
        
        /* Garantir que o ícone está centrado */
        .sidebar-toggle i {
          font-size: 16px !important;
        }
      "))
    )
    ,
    
    
    tabItems(
      # Tab principal de análise
      tabItem(tabName = "main_analysis",
              
              fluidRow(
                # Resumo estatístico
                box(
                  title = "Resumo da Análise", 
                  status = "primary", 
                  solidHeader = TRUE,
                  width = 12,
                  collapsible = TRUE,
                  
                  fluidRow(
                    column(3, 
                           valueBoxOutput("total_genes", width = 12)
                    ),
                    column(3, 
                           valueBoxOutput("de_genes", width = 12)
                    ),
                    column(3, 
                           valueBoxOutput("upregulated", width = 12)
                    ),
                    column(3, 
                           valueBoxOutput("downregulated", width = 12)
                    )
                  ),
                  
                  br(),
                  verbatimTextOutput("analysis_summary")
                )
              ),
              
              fluidRow(
                # Volcano Plot - agora ocupa mais espaço
                box(
                  title = "Volcano Plot", 
                  status = "primary", 
                  solidHeader = TRUE,
                  width = 8,
                  collapsible = TRUE,
                  
                  plotlyOutput("volcano_plot", height = "600px")
                ),
                
                # PCA Plot
                box(
                  title = "Análise PCA", 
                  status = "primary", 
                  solidHeader = TRUE,
                  width = 4,
                  collapsible = TRUE,
                  
                  plotlyOutput("pca_plot", height = "600px")
                )
              ),
              
              fluidRow(
                # Tabela de genes DE - agora tem mais destaque
                box(
                  title = "Top 20 Genes Diferencialmente Expressos", 
                  status = "info", 
                  solidHeader = TRUE,
                  width = 12,
                  collapsible = TRUE,
                  
                  div(style = "height: 400px; overflow-y: auto; overflow-x: auto; width: 100%;",
                      tableOutput("de_genes_table")
                  )
                )
              ),
              
              fluidRow(
                # Boxplots exploratórios - mantidos
                box(
                  title = "Distribuição das Contagens (Boxplots)", 
                  status = "primary", 
                  solidHeader = TRUE,
                  width = 12,
                  collapsible = TRUE,
                  
                  plotOutput("boxplot_analysis", height = "1500px")
                )
              )
      ),
      # Nova aba de heatmaps
      tabItem(tabName = "heatmaps",
              
              fluidRow(
                box(
                  title = "Heatmap Dinâmico - Comparação Ativa", 
                  status = "primary", 
                  solidHeader = TRUE,
                  width = 12,
                  collapsible = TRUE,
                  
                  p("Este heatmap apresenta os genes mais diferencialmente expressos na comparação selecionada, permitindo observar, de forma global, como estes genes se comportam em todas as amostras e condições do estudo."),
                  
                  plotOutput("heatmap_dynamic", height = "700px")
                )
              ),
              
              
              fluidRow(
                box(
                  title = "Heatmap Global - DMSO vs MW por Time-point", 
                  status = "warning", 
                  solidHeader = TRUE,
                  width = 12,
                  collapsible = TRUE,
                  
                  p("Este heatmap apresenta os genes mais diferencialmente expressos entre todas as amostras,permitindo observar, de forma global, como estes genes se comportam em todas as amostras e condições do estudo."),
                  
                  plotOutput("heatmap_comparison", height = "800px")
                )
              )
      ),
      # Tab de  informação
      tabItem(tabName = "methods",
              fluidRow(
                box(
                  title = "Informações sobre a APP", 
                  status = "info", 
                  solidHeader = TRUE,
                  width = 12,
                  
                  p("Esta aplicação analisa dados de RNA-seq de células MCF7 (linha celular de cancro da mama) 
              tratadas com uma combinação de Mefloquina e Withaferina A (MW) em três time-points 
              (0h, 9h, 24h), comparadas com controlo DMSO."),
                  
                  h3("Parâmetros"),
                  
                  # Título simples
                  tags$p(tags$b("Escolher Comparação:")),
                  
                  # Bullets com explicações
                  tags$ul(
                    tags$li("MW_24h vs MW_0h: Para estudar o efeito do tempo no tratamento."),
                    tags$li("MW_24h vs DMSO_24h: Para observar o efeito do tratamento no mesmo time-point."),
                    tags$li("MW_9h vs DMSO_9h: Análise intermédia."),
                    tags$li("MW_0h vs DMSO_0h: Controlo basal, útil para verificar se os grupos são comparáveis antes do tratamento.")
                  ),
                  
                  # Título simples
                  tags$p(tags$b("Thresholds:")),
                  
                  # Bullets com explicações
                  tags$ul(
                    tags$li("Log2FoldChange > 0: Gene sobre-expresso na condição teste."),
                    tags$li("Log2FoldChange < 0: Gene sub-expresso na condição teste."),
                    tags$li("padj < 0.05: Gene significativamente diferencial (após correção FDR).")
                  ),
                  
                  # Título simples
                  tags$p(tags$b("Número de genes no heatmap:")),
                  
                  # Bullet com explicação
                  tags$ul(
                    tags$li("Quantidade de genes diferencialmente expressos a visualizar no heatmap.")
                  )
                  ,
                  
                  h3("Pipeline de Análise"),
                  tags$ol(
                    tags$li("Carregar e processar dados de contagem"),
                    tags$li("Filtrar genes com baixa expressão"),
                    tags$li("Normalizar com DESeq2 (median ratio method)"),
                    tags$li("Estimativa de dispersão"),
                    tags$li("Teste estatístico (Wald test)"),
                    tags$li("Correção para múltiplas comparações (FDR)")
                  ),
                  
                  h3("Ficheiros Necessários"),
                  tags$ul(
                    tags$li("gene_count_matrix.csv ou transcript_count_matrix.csv"),
                    tags$li("sample_info.txt"),
                    tags$li("transcripts_annot.txt (se usar transcript_count_matrix)")
                  ),
                  
                  br(),
                  p(em("Pipeline baseado no trabalho do Prof. Jorge Cabral, adaptado de Andreia Reis, 
                 Margarida Ferreira e Rita Guimarães"))
                )
              )
      )
    )
  )
)

# ============================================================================
# LÓGICA DO SERVIDOR
# ============================================================================
options(shiny.error = function() {
  traceback(20)
})



server <- function(input, output, session) {
  
  # Dados reativos
  data_processed <- reactive({
    withProgress(message = 'Carregando dados...', {
      incProgress(0.3, detail = "Preparando dataset...")
      prepare_dataset()
    })
  })
  
  # Análise DESeq2 reativa
  deseq_results <- reactive({
    
    req(data_processed(), input$comparison_type)
    
    withProgress(message = 'Executando análise DESeq2...', {
      
      data <- data_processed()
      
      # Determinar referência e experimental baseado na seleção
      comp_map <- list(
        "MW_24_vs_MW_0" = list(ref = "MW_0", exp = "MW_24"),
        "MW_24_vs_DMSO_24" = list(ref = "DMSO_24", exp = "MW_24"),
        "MW_9_vs_DMSO_9" = list(ref = "DMSO_9", exp = "MW_9"),
        "MW_0_vs_DMSO_0" = list(ref = "DMSO_0", exp = "MW_0")
      )
      
      comparison <- comp_map[[input$comparison_type]]
      
      incProgress(0.7, detail = "Executando DESeq2...")
      
      run_deseq2_analysis(
        data$cts, 
        data$coldata, 
        data$variable_of_interest,
        comparison$ref,
        comparison$exp
      )
    })
  })
  
  # Preparar dados dos resultados
  results_data <- reactive({
    req(deseq_results())
    
    res <- deseq_results()$res
    
    # Verificar se temos resultados válidos
    if(is.null(res) || nrow(res) == 0) {
      return(NULL)
    }
    
    res_df <- as.data.frame(res)
    res_df$gene_id <- rownames(res_df)
    
    # Debug
    cat("Processando resultados DESeq2:\n")
    cat("- Número de genes:", nrow(res_df), "\n")
    cat("- Primeiros gene_ids:", paste(head(res_df$gene_id, 3), collapse = ", "), "\n")
    
    # Extrair nome do gene se houver anotação
    if(any(grepl("\\|", res_df$gene_id))) {
      res_df$gene_name <- sapply(strsplit(res_df$gene_id, "\\|"), function(x) {
        if(length(x) > 1) x[2] else x[1]
      })
    } else {
      res_df$gene_name <- res_df$gene_id
    }
    
    # Remover linhas com NA em padj para evitar problemas
    res_df <- res_df[!is.na(res_df$padj), ]
    
    # Classificar significância
    res_df$significant <- ifelse(res_df$padj < input$pval_threshold & 
                                   abs(res_df$log2FoldChange) > input$lfc_threshold,
                                 "Significativo", "Não significativo")
    
    res_df$regulation <- ifelse(res_df$log2FoldChange > 0, "Sobre-expresso", "Sub-expresso")
    
    cat("- Genes significativos:", sum(res_df$significant == "Significativo"), "\n")
    
    res_df
  })
  
  # ============================================================================
  # OUTPUTS REATIVOS
  # ============================================================================
  
  # Value boxes
  output$total_genes <- renderValueBox({
    valueBox(
      value = ifelse(is.null(results_data()), 0, nrow(results_data())),
      subtitle = "Total de Genes Analisados",
      icon = icon("dna"),
      color = "blue"
    )
  })
  
  output$de_genes <- renderValueBox({
    valueBox(
      value = ifelse(is.null(results_data()), 0, 
                     sum(results_data()$padj < input$pval_threshold, na.rm = TRUE)),
      subtitle = "Genes Diferencialmente Expressos",
      icon = icon("chart-line"),
      color = "green"
    )
  })
  
  output$upregulated <- renderValueBox({
    valueBox(
      value = ifelse(is.null(results_data()), 0,
                     sum(results_data()$padj < input$pval_threshold & 
                           results_data()$log2FoldChange > input$lfc_threshold, na.rm = TRUE)),
      subtitle = "Genes Sobre-expressos",
      icon = icon("arrow-up"),
      color = "red"
    )
  })
  
  output$downregulated <- renderValueBox({
    valueBox(
      value = ifelse(is.null(results_data()), 0,
                     sum(results_data()$padj < input$pval_threshold & 
                           results_data()$log2FoldChange < -input$lfc_threshold, na.rm = TRUE)),
      subtitle = "Genes Sub-expressos",
      icon = icon("arrow-down"),
      color = "yellow"
    )
  })
  
  # Resumo da análise
  output$analysis_summary <- renderText({
    req(deseq_results())
    
    comparison_labels <- list(
      "MW_24_vs_MW_0" = "MW 24h vs MW 0h (Efeito temporal do tratamento)",
      "MW_24_vs_DMSO_24" = "MW 24h vs DMSO 24h (Efeito do tratamento às 24h)",
      "MW_9_vs_DMSO_9" = "MW 9h vs DMSO 9h (Efeito do tratamento às 9h)",
      "MW_0_vs_DMSO_0" = "MW 0h vs DMSO 0h (Controlo basal)"
    )
    
    paste0(
      "Comparação Ativa: ", comparison_labels[[input$comparison_type]], "\n",
      "Threshold p-value ajustado: ", input$pval_threshold, "\n",
      "Threshold Log2FoldChange: ±", input$lfc_threshold, "\n",
      "Método de correção: FDR (False Discovery Rate)\n",
      "Nome da comparação: ", deseq_results()$comparison_name
    )
  })
  
  # Volcano Plot
  output$volcano_plot <- renderPlotly({
    req(results_data())
    
    res_df <- results_data()
    res_df <- res_df[!is.na(res_df$padj), ]
    
    # Top 10 genes para anotação
    top_genes <- res_df[order(res_df$padj), ][1:10, ]
    
    p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(aes(color = significant, text = paste("Gene:", gene_name,
                                                       "<br>Log2FC:", round(log2FoldChange, 3),
                                                       "<br>p-adj:", format(padj, scientific = TRUE))),
                 alpha = 0.7, size = 1.5) +
      scale_color_manual(values = c("Significativo" = "#e74c3c", "Não significativo" = "#95a5a6")) +
      geom_vline(xintercept = c(-input$lfc_threshold, input$lfc_threshold), 
                 linetype = "dashed", color = "gray50") +
      geom_hline(yintercept = -log10(input$pval_threshold), 
                 linetype = "dashed", color = "gray50") +
      theme_minimal() +
      labs(title = "Volcano Plot - Expressão Génica Diferencial",
           x = "Log2 Fold Change",
           y = "-Log10(p-value ajustado)",
           color = "Significância") +
      theme(legend.position = "bottom")
    
    # Adicionar labels para top genes
    if(nrow(top_genes) > 0) {
      p <- p + geom_text_repel(data = top_genes,
                               aes(label = gene_name),
                               size = 3,
                               max.overlaps = 10)
    }
    
    ggplotly(p, tooltip = "text") %>%
      layout(legend = list(orientation = "h", x = 0.3, y = -0.2))
  })
  
  # PCA Plot
  output$pca_plot <- renderPlotly({
    req(deseq_results())
    
    vsd <- deseq_results()$vsd
    data <- data_processed()
    
    # Preparar dados para PCA
    pca_data <- plotPCA(vsd, intgroup = colnames(data$coldata)[3], returnData = TRUE)
    
    # Adicionar informações de tratamento e tempo
    pca_data$treatment <- sapply(strsplit(as.character(pca_data$group), "_"), `[`, 1)
    pca_data$timepoint <- sapply(strsplit(as.character(pca_data$group), "_"), `[`, 2)
    
    p <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
      geom_point(aes(color = treatment, shape = timepoint, 
                     text = paste("Amostra:", name,
                                  "<br>Tratamento:", treatment,
                                  "<br>Tempo:", timepoint, "h")),
                 size = 4, alpha = 0.8) +
      scale_color_brewer(type = "qual", palette = "Set1") +
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      theme_minimal() +
      labs(title = "Análise de Componentes Principais (PCA)",
           color = "Tratamento",
           shape = "Tempo (h)") +
      theme(legend.position = "bottom")
    
    ggplotly(p, tooltip = "text") %>%
      layout(legend = list(orientation = "h", x = 0.2, y = -0.2))
  })
  
  # Tabela de genes DE
  output$de_genes_table <- renderTable({
    req(results_data())
    
    res_df <- results_data()
    sig_genes <- res_df[!is.na(res_df$padj) & res_df$padj < input$pval_threshold, ]
    sig_genes <- sig_genes[order(sig_genes$padj), ]
    
    if(nrow(sig_genes) > 0) {
      top_20 <- head(sig_genes, 20)
      
      table_data <- data.frame(
        "Gene" = top_20$gene_name,
        "Log2FC" = round(top_20$log2FoldChange, 3),
        "p-value" = format.pval(top_20$pvalue, digits = 3, eps = .Machine$double.eps),
        "p-adj" = format.pval(top_20$padj, digits = 3, eps = .Machine$double.eps),
        "Regulação" = ifelse(top_20$log2FoldChange > 0, "↑", "↓")
      )
      
      # Transpor a tabela para ter genes nas colunas
      t_table <- t(table_data)
      
      # Ajustar nomes das colunas para os genes
      colnames(t_table) <- t_table[1, ] # primeira linha é gene
      t_table <- t_table[-1, ] # remover linha dos nomes
      
      # Converter em data.frame e manter os nomes das linhas
      t_table_df <- as.data.frame(t_table, stringsAsFactors = FALSE)
      
      # Alterar o nome da coluna da variável para "Gene"
      t_table_df <- cbind(Gene = rownames(t_table_df), t_table_df)
      rownames(t_table_df) <- NULL
      
      t_table_df
    } else {
      data.frame("Mensagem" = "Nenhum gene significativo encontrado")
    }
  }, striped = TRUE, hover = TRUE, spacing = "xs")
  
  # Boxplots
  output$boxplot_analysis <- renderPlot({
    req(deseq_results())
    
    dds <- deseq_results()$dds
    vsd <- deseq_results()$vsd
    data <- data_processed()
    
    # Preparar dados para boxplots
    plot_data <- list(
      "Contagens Brutas (log10)" = log10(as.matrix(counts(dds)) + 1),
      "Contagens Normalizadas (log10)" = log10(as.matrix(counts(dds, normalized = TRUE)) + 1),
      "VST (Variance Stabilizing Transformation)" = as.matrix(assay(vsd))
    )
    
    # Cores por tratamento
    treatment_colors <- c("DMSO_0" = "#f1c40f", "DMSO_9" = "#f39c12", "DMSO_24" = "#e67e22",
                          "MW_0" = "#3498db", "MW_9" = "#2980b9", "MW_24" = "#1abc9c")
    
    col_colors <- treatment_colors[data$coldata[, 3]]
    
    # Define o layout uma única vez para os 3 plots em coluna
    par(mfrow = c(3, 1))
    
    for(i in 1:length(plot_data)) {
      num_boxes <- ncol(plot_data[[i]])
      positions <- seq(1, num_boxes * 3, by = 3)
      
      # Ajusta margens para cada gráfico para evitar sobreposição vertical
      if(i == 1) {
        par(mar = c(10, 5, 5, 7))  # gráfico 1: margem superior maior para legenda
      } else if(i == 2) {
        par(mar = c(10, 5, 5, 7))   # gráfico 2: margens médias
      } else {
        par(mar = c(10, 5, 5, 7))   # gráfico 3: margem superior menor
      }
      
      boxplot(plot_data[[i]],
              at = positions,
              xlim = c(0, max(positions) + 3),
              main = names(plot_data)[i],
              xaxt = "n",
              las = 2,
              col = col_colors,
              outline = FALSE,
              cex.axis = 1.5,
              cex.lab = 1.5,
              cex.main = 1.8,
              ylab = "Expressão")
      
      axis(1, at = positions, labels = colnames(plot_data[[i]]), las = 2, cex.axis = 1.5)
      
      legend("topright", 
             legend = names(treatment_colors),
             fill = treatment_colors,
             cex = 1.2,
             bty = "n",
             title = "Condições")
    }
  })
  # Heatmap dinâmico (baseado na comparação ativa)
  output$heatmap_dynamic <- renderPlot({
    req(results_data(), deseq_results())
    
    # Obter top genes da comparação atual
    res_df <- results_data()
    sig_genes <- res_df[!is.na(res_df$padj) & res_df$padj < input$pval_threshold, ]
    
    if(nrow(sig_genes) == 0) {
      plot.new()
      text(0.5, 0.5, "Nenhum gene significativo encontrado\nna comparação atual", 
           cex = 1.5, col = "red")
      return()
    }
    
    # Top genes por p-valor ajustado
    top_genes <- sig_genes[order(sig_genes$padj), ][1:min(input$n_genes_heatmap, nrow(sig_genes)), ]$gene_id
    
    # Dados transformados (VST)
    vsd <- deseq_results()$vsd
    
    # Verificar se os genes existem nos dados
    available_genes <- intersect(top_genes, rownames(vsd))
    
    if(length(available_genes) == 0) {
      plot.new()
      text(0.5, 0.5, "Genes selecionados não encontrados nos dados", 
           cex = 1.5, col = "red")
      return()
    }
    
    mat <- assay(vsd)[available_genes, ]
    
    # Preparar nomes das colunas
    data <- data_processed()
    col_annotations <- data.frame(
      Tratamento = sapply(strsplit(as.character(data$coldata[, 3]), "_"), `[`, 1),
      Tempo = paste0(sapply(strsplit(as.character(data$coldata[, 3]), "_"), `[`, 2), "h")
    )
    rownames(col_annotations) <- colnames(mat)
    
    # Cores para anotações
    ann_colors <- list(
      Tratamento = c("DMSO" = "#3498db", "MW" = "#e74c3c"),
      Tempo = c("0h" = "#f1c40f", "9h" = "#f39c12", "24h" = "#e67e22")
    )
    
    # Criar heatmap
    pheatmap(mat, 
             scale = "row",
             clustering_distance_rows = "correlation",
             clustering_distance_cols = "euclidean",
             annotation_col = col_annotations,
             annotation_colors = ann_colors,
             main = paste("Top", length(available_genes), "genes -", input$comparison_type),
             fontsize = 10,
             fontsize_row = 8,
             fontsize_col = 8,
             show_rownames = length(available_genes) <= 30,
             show_colnames = TRUE)
  })
  
  
  
  # Heatmap comparativo geral (manter igual)
  output$heatmap_comparison <- renderPlot({
    req(data_processed(), deseq_results())
    
    data <- data_processed()
    vsd <- deseq_results()$vsd
    
    # Usar genes mais variáveis para comparação geral
    gene_vars <- apply(assay(vsd), 1, var)
    gene_vars <- gene_vars[!is.na(gene_vars) & gene_vars > 0]
    var_genes <- names(sort(gene_vars, decreasing = TRUE)[1:min(input$n_genes_heatmap, length(gene_vars))])
    
    mat <- assay(vsd)[var_genes, ]
    
    # Anotações completas
    col_annotations <- data.frame(
      Tratamento = sapply(strsplit(as.character(data$coldata[, 3]), "_"), `[`, 1),
      Tempo = paste0(sapply(strsplit(as.character(data$coldata[, 3]), "_"), `[`, 2), "h")
    )
    rownames(col_annotations) <- colnames(mat)
    
    # Cores para anotações
    ann_colors <- list(
      Tratamento = c("DMSO" = "#3498db", "MW" = "#e74c3c"),
      Tempo = c("0h" = "#f1c40f", "9h" = "#f39c12", "24h" = "#e67e22")
    )
    
    pheatmap(mat, 
             scale = "row",
             clustering_distance_rows = "correlation",
             clustering_distance_cols = "euclidean",
             annotation_col = col_annotations,
             annotation_colors = ann_colors,
             main = "Comparação Geral - DMSO vs MW por Time-point",
             fontsize = 9,
             fontsize_row = 7,
             show_rownames = length(var_genes) <= 30,
             show_colnames = TRUE)
  })
  # Download handler
  output$download_results <- downloadHandler(
    filename = function() {
      paste0("resultados_deseq2_", input$comparison_type, "_", 
             format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      req(results_data())
      
      res_df <- results_data()
      sig_genes <- res_df[!is.na(res_df$padj) & res_df$padj < input$pval_threshold, ]
      sig_genes <- sig_genes[order(sig_genes$padj), ]
      
      export_data <- data.frame(
        Gene_ID = sig_genes$gene_id,
        Gene_Name = sig_genes$gene_name,
        Log2FoldChange = sig_genes$log2FoldChange,
        pvalue = sig_genes$pvalue,
        padj = sig_genes$padj,
        Regulation = sig_genes$regulation,
        baseMean = sig_genes$baseMean,
        lfcSE = sig_genes$lfcSE,
        stat = sig_genes$stat
      )
      
      write.csv(export_data, file, row.names = FALSE)
    }
  )
}  
#============================================================================
#FIM DO SERVIDOR
#============================================================================
shinyApp(ui = ui, server = server)