library(shiny)
library(bslib)
library(shinythemes)
library(scales)  # Para manipulação de cores
library(dplyr)
library(stringr)
library(tools)  # Para manipulação de arquivos
library(rmarkdown)
library(knitr)

#bslib/shinythemes: para temas e estilização
#scales/dplyr/stringr: para manipulação de dados
#shiny: framework principal
#rmarkdown/knitr: para geração de relatórios


# Gera relatório PDF com dados do paciente
# Usa template RMarkdown
# Inclui tratamento de erros
generate_pdf_report <- function(input_data) {
  cat("Input Data for PDF Generation:\n")
  print(str(input_data))
  
  # Função interna para validar e converter dados
  check_and_convert <- function(value, default = "Não Informado", type = "character") {
    if (is.null(value) || is.na(value)) {
      return(default)
    }
    #converter o valor para o tipo especificado
    tryCatch({
      switch(type,
             "character" = as.character(value),
             "numeric" = format(as.numeric(value), nsmall = 0),
             "date" = format(as.Date(value), "%d de %B de %Y"),
             as.character(value))
    }, error = function(e) {
      #Exibe mensagem de erro em caso de falha na conversão
      cat("Erro na conversão:", conditionMessage(e), "\n")
      return(default)
    })
  }
  
  
  # Prepara os dados do paciente com validações
  prepared_data <- list(
    date = check_and_convert(input_data$date %||% Sys.Date(), 
                             default = format(Sys.Date(), "%d de %B de %Y"), 
                             type = "date"),
    name = check_and_convert(input_data$name),
    age = check_and_convert(input_data$age, default = "0", type = "numeric"),
    cc_number = check_and_convert(input_data$cc_number),
    utente_number = check_and_convert(input_data$utente_number),
    fev1 = check_and_convert(input_data$fev1, default = "0", type = "numeric"),
    gold_classification = check_and_convert(input_data$gold_classification),
    exacerbations_hosp = check_and_convert(input_data$exacerbations_hosp, default = "0", type = "numeric"),
    exacerbations_no_hosp = check_and_convert(input_data$exacerbations_no_hosp, default = "0", type = "numeric"),
    symptom_details = check_and_convert(input_data$symptom_details),
    treatment_suggestion = check_and_convert(input_data$treatment_suggestion)
  )
  
  # Criar arquivo RMarkdown temporário
  temp_rmd <- tempfile(fileext = ".Rmd")
  temp_pdf <- tempfile(fileext = ".pdf")
  # Define o conteúdo do relatório usando o template RMarkdown
  report_content <- sprintf("
---
title: 'Relatório de Avaliação DPOC'
date: '%s'
output: 
  pdf_document:
    latex_engine: xelatex
---

## Dados do Paciente

- **Nome:** %s
- **Idade:** %s anos
- **Cartão de Cidadão:** %s
- **Número de Utente:** %s
- **Data da Avaliação:** %s

## Avaliação Clínica DPOC

### Classificação de Obstrução
- **FEV1:** %s%%
- **Classificação GOLD:** %s

### Histórico de Exacerbações
- **Exacerbações com Hospitalização:** %s
- **Exacerbações sem Hospitalização:** %s

### Questionário de Sintomas
%s

### Recomendação de Tratamento
%s
", 
                            prepared_data$date,
                            prepared_data$name,
                            prepared_data$age,
                            prepared_data$cc_number,
                            prepared_data$utente_number, 
                            prepared_data$date,
                            prepared_data$fev1,
                            prepared_data$gold_classification,
                            prepared_data$exacerbations_hosp,
                            prepared_data$exacerbations_no_hosp,
                            prepared_data$symptom_details,
                            prepared_data$treatment_suggestion)
  
  # Tratamento de erros durante a escrita e renderização
  tryCatch({
    # Escreve conteúdo no arquivo temporário
    writeLines(report_content, temp_rmd)
    
    # Renderiza o RMarkdown como PDF
    rmarkdown::render(
      temp_rmd, 
      output_file = temp_pdf, 
      output_format = "pdf_document",
      clean = TRUE,  # Limpa arquivos intermediários
      quiet = FALSE  # Mostra mensagens de erro
    )
    
    # Verificar se o arquivo PDF foi criado
    if (!file.exists(temp_pdf)) {
      stop("Falha na geração do PDF")
    }
    # Retorna o caminho do arquivo PDF
    return(temp_pdf)
  }, error = function(e) {
    # Log detalhado do erro
    cat("Erro na geração do PDF: ", conditionMessage(e), "\n")
    cat("Detalhes do erro:\n")
    print(e)
    # Interrompe a execução e informa o erro
    stop(paste("Erro na geração do relatório:", conditionMessage(e)))
  }, finally = {
    # Remove o arquivo temporário RMarkdown, independentemente 
    unlink(temp_rmd)
  })
}



# Função para validar e limpar entradas
validate_input <- function(input, type = "text", max_length = NULL) {
  # Verifica se o input está vazio ou é nulo. Caso positivo, retorna NULL
  if (is.null(input) || input == "") return(NULL)
  # Usa o switch para diferentes tipos de validação com base no tipo do dado
  switch(type,
         "text" = {
           # Remove caracteres especiais, mantém letras, espaços e acentos
           cleaned <- str_replace_all(input, "[^[:alnum:][:space:]áéíóúãõâêîôûàèìòùäëïöüçÁÉÍÓÚÃÕÂÊÎÔÛÀÈÌÒÙÄËÏÖÜÇ]", "") %>%
             str_trim() # Remove espaços em branco no início e no final
           
           # Limita o comprimento do texto se 'max_length' for especificado
           if (!is.null(max_length)) cleaned <- substr(cleaned, 1, max_length)
           return(cleaned) # Retorna o texto limpo
         },
         "numeric" = {
           # Converte para numérico, remove não-números
           cleaned <- as.numeric(str_replace_all(input, "[^0-9.]", ""))
           return(ifelse(is.na(cleaned), 0, cleaned))
           # Se a conversão para numérico falhar (NA), retorna 0
         },
         "cc_number" = {
           # Remove tudo que não for número ou letra maiúscula e limita a 12 caracteres
           cleaned <- str_replace_all(input, "[^A-Z0-9]", "") %>%
             substr(1, 12)  # Mantém apenas os primeiros 12 caracteres
           return(cleaned)  # Retorna o valor limpo
         },
         "utente_number" = {
           # Remove tudo exceto números, limita ao comprimento do número de utente
           cleaned <- str_replace_all(input, "[^0-9]", "") %>%
             substr(1, 9) # Mantém apenas os primeiros 9 caracteres
           return(cleaned)
         }
  )
}
# Função para gerar cores baseadas na classificação GOLD
get_gold_color <- function(fev1) {
  if (is.null(fev1) || fev1 > 100) return("#CCCCCC")  # Cinza para valores inválidos
  
  # Escala de cores do verde (menos grave) ao vermelho (mais grave)
  if (fev1 >= 80) {
    return("#2ecc71")  # Verde (GOLD 1 - Leve)
  } else if (fev1 >= 50) {
    return("#f1c40f")  # Amarelo (GOLD 2 - Moderado)
  } else if (fev1 >= 30) {
    return("#e67e22")  # Laranja (GOLD 3 - Grave)
  } else {
    return("#e74c3c")  # Vermelho (GOLD 4 - Muito Grave)
  }
}


#Cria uma interface com barra de navegação
#Inclui logo e título
#Usa tema personalizado com cores GOLD (dourado)
#Dividido em abas: Introdução, Dados do Paciente, Classificação da Obstrução, etc.
ui <- navbarPage(
  #theme = shinytheme("cerulean"),
  #title = "Ferramenta de Avaliação GOLD ABE",
  title = div(
    tags$img(src = "gold.png", height = "50px", style = "margin-right: 10px;"),
    "Ferramenta de Avaliação GOLD ABE"
  ),
  id = "main_tabset",
  
  theme = bs_theme(
    version = 3, 
    primary = "#e8c251", 
    #navbar_bg = "#FFD700", 
    #navbar_color = "#000000"
    heading_font = bslib::font_google("Montserrat"), 
    base_font = bslib::font_google("Roboto") 
  ),
  
  tags$head(
    tags$style(HTML("
      .navbar-default { background-color: #e8c251 !important; border-color: #e8c251 !important; }
      .navbar-default .navbar-brand { color: #000000 !important; font-weight: bold; }
      .navbar-default .navbar-nav > li > a { color: #000000 !important; }
      .navbar-default .navbar-nav > li > a:hover { background-color: #fff !important; }
      h1, h2, h3 { color: #e8c251; font-weight: bold; }
      .btn-primary { background-color: #e8c251; color: #000000; border: none; }
    "))
  ),
  
  # Aba 1: Introdução
  tabPanel(
    "Introdução",
    h1("Aplicação GOLD 2024"),
    p("Esta aplicação foi desenvolvida para ajudar os clínicos na escolha da 
         abordagem farmacológica inicial para pacientes com DPOC."),
    p("Autores: grupo B"),
    p(em("Nota: Esta aplicação é apenas uma ferramenta auxiliar e não substitui 
            o julgamento clínico. Os autores não se responsabilizam por decisões 
            baseadas exclusivamente nesta aplicação."))
  ),
  
  # Aba 2: Dados do Paciente
  tabPanel(
    "Dados do Paciente",
    fluidPage(
      h3("Dados do Paciente"),
      textInput("name", "Nome do Paciente:", value = ""),
      numericInput("age", "Idade do Paciente:", value = 50, min = 18, max = 120),
      textInput("address", "Morada do Paciente:", value = ""),
      textInput("cc_number", "Número do Cartão de Cidadão:", value = ""),
      textInput("utente_number", "Número de Utente:", value = ""),
      dateInput("date", "Data do Questionário:", value = Sys.Date())
    )
  ),
  
  # Aba 3: Classificação da Obstrução
  tabPanel(
    "Classificação da Obstrução",
    h2("Classificação da Severidade da Obstrução (Pós-broncodilatador)"),
    numericInput("fev1", "Insira o valor de FEV1 (% do previsto):", 
                 value = 0, min = 0, max = 100, step = 0.1),
    # Novo elemento para visualização colorida
    uiOutput("gold_classification_visual")
  ),
  
  
  # Aba 4: Tratamento Farmacológico
  tabPanel(
    "Tratamento Farmacológico",
    fluidPage(
      # Linha principal: Histórico de Exacerbações e Outros Fatores
      fluidRow(
        # Histórico de Exacerbações
        column(
          width = 6,
          div(
            h3("Histórico de Exacerbações", style = "text-align: left; margin-left: 400px;"),
            div(
              numericInput(
                "exacerbations_hosp", 
                "Número de exacerbações COM hospitalização (último ano):", 
                value = NA, min = 0, max = 100
              ),
              style = "margin-left: 400px; margin-bottom: 10px;"
            ),
            div(
              numericInput(
                "exacerbations_no_hosp", 
                "Número de exacerbações SEM hospitalização (último ano):", 
                value = NA, min = 0, max = 100
              ),
              style = "margin-left: 400px; margin-bottom: 20px;"
            )
          )
        ),
        # Outros Fatores (Condicional)
        column(
          width = 6,
          conditionalPanel(
            condition = "output.show_eosinophils_asthma",
            div(
              h3("Outros Fatores", style = "text-align: left; margin-left: 20px;"),
              div(
                numericInput(
                  "blood_eosinophils", 
                  "Contagem de eosinófilos sanguíneos (células/μL):", 
                  value = NA, min = 0, max = 1000
                ),
                style = "margin-left: 20px; margin-bottom: 10px;"
              ),
              div(
                checkboxInput("asthma", "Histórico de Asma", value = FALSE),
                style = "margin-left: 20px;"
              )
            )
          )
        )
      ),
      # Linha de separação
      hr(style = "border: 1px solid #ddd; margin-top: 20px; margin-bottom: 20px;"),
      # Linha inferior: Questionário de Sintomas
      fluidRow(
        column(
          width = 12,
          div(
            h3("Questionário de Sintomas", style = "text-align: left; margin-left: 600px;"),
            div(
              radioButtons(
                "questionnaire_choice", 
                "Escolha o questionário:", 
                choices = c("mMRC" = "mMRC", "CAT" = "CAT")
              ),
              style = "margin-left: 600px; margin-top: 10px;"
            ),
            div(
              uiOutput("questionnaire_ui"),
              style = "margin-left: 600px; margin-top: 10px;"
            )
          )
        )
      )
    )
  ),
  
  
  
  # Aba 4: Tratamento Farmacológico
  tabPanel(
    "Resultados e Relatório",
    fluidPage(
      # Campo para os resultados que aparece dinamicamente
      fluidRow(
        column(
          width = 12,
          h3("Resultados e Relatório"),
          verbatimTextOutput("summary"),
          verbatimTextOutput("gold_classification_resumo"),
          verbatimTextOutput("group"),
          verbatimTextOutput("suggestion"),
          downloadButton("download_pdf", "Download Relatório em PDF", class = "btn-primary")
        )
      )
    )
  )
)
#Contém toda a lógica reativa da aplicação
#Processa inputs do usuário
#Atualiza a interface em tempo real
#Gerencia cálculos e validações
server <- function(input, output, session) {
  gold_grade <- reactiveVal(NA) # Variável reativa para armazenar a classificação GOLD
  
  # Renderiza dinamicamente o questionário escolhido pelo clínico (CAT ou mMRC)
  output$questionnaire_ui <- renderUI({
    if (input$questionnaire_choice == "CAT") {
      fluidPage(
        sliderInput("q1", "Eu nunca tusso ↔ Eu tusso o tempo todo", min = 0, max = 5, value = 0),
        sliderInput("q2", "Eu não tenho muco no peito ↔ O meu peito está completamente cheio de muco", min = 0, max = 5, value = 0),
        sliderInput("q3", "O meu peito não se sente apertado ↔ O meu peito sente-se muito apertado", min = 0, max = 5, value = 0),
        sliderInput("q4", "Eu não fico sem ar ao subir um lançe de escadas ↔ Eu fico sem ar ao subir um lançe de escadas", min = 0, max = 5, value = 0),
        sliderInput("q5", "Eu não sou limitado a fazer tarefas domésticas ↔ Eu sou muito limitado a fazer tarefas domésticas", min = 0, max = 5, value = 0),
        sliderInput("q6", "Eu estou confiante ao sair de casa ↔ Eu não estou confiante ao sair de casa devido à minha condição pulmonar", min = 0, max = 5, value = 0),
        sliderInput("q7", "Eu durmo profundamente ↔ Não durmo profundamente devido ao meu problema pulmonar", min = 0, max = 5, value = 0),
        sliderInput("q8", "Eu tenho muita energia ↔ Não tenho energia alguma", min = 0, max = 5, value = 0)
      )
    } else if (input$questionnaire_choice == "mMRC") {
      radioButtons("mmrc", "Escolha a descrição que melhor se aplica:", 
                   choices = list(
                     "Sem dispneia" = 0,
                     "Falta de ar ao caminhar rápido ou subir leve declive" = 1,
                     "Precisa parar para respirar ao caminhar em terreno plano" = 2,
                     "Falta de ar após 100 metros ou poucos minutos" = 3,
                     "Falta de ar ao vestir-se ou em atividades diárias" = 4
                   ))
    }
  })
  ### Validação e atualização dos inputs do usuário
  ### Classificação da Obstrução
  observe({
    # Input validations
    updateTextInput(session, "name", 
                    value = validate_input(input$name, type = "text", max_length = 100))
    
    updateTextInput(session, "address", 
                    value = validate_input(input$address, type = "text", max_length = 200))
    
    updateTextInput(session, "cc_number", 
                    value = validate_input(input$cc_number, type = "cc_number"))
    
    updateTextInput(session, "utente_number", 
                    value = validate_input(input$utente_number, type = "utente_number"))
    
    updateNumericInput(session, "age", 
                       value = max(18, min(input$age, 120)))
    
    updateNumericInput(session, "fev1", 
                       value = max(0, min(input$fev1, 100)))
    
    updateNumericInput(session, "exacerbations_hosp", 
                       value = max(0, min(round(input$exacerbations_hosp), 100)))
    
    updateNumericInput(session, "exacerbations_no_hosp", 
                       value = max(0, min(round(input$exacerbations_no_hosp), 100)))
    
    updateNumericInput(session, "blood_eosinophils", 
                       value = max(0, min(input$blood_eosinophils, 1000)))
    
    # GOLD Grade Calculation
    current_gold_grade <- ifelse(is.null(input$fev1) || input$fev1 > 100, NA,
                                 ifelse(input$fev1 >= 80, "GOLD 1 - Leve ",
                                               ifelse(input$fev1 >= 50, "GOLD 2 - Moderado",
                                                      ifelse(input$fev1 >= 30, "GOLD 3 - Grave", "GOLD 4 - Muito Grave"))))
    
    
    
    
    
    gold_grade(current_gold_grade)  # Atualiza o valor reativo
    # Atualiza a visualização da classificação GOLD com base na cor
    output$gold_classification_visual <- renderUI({
      
      req(input$fev1)# Verifica se FEV1 foi inserido
      
      color <- get_gold_color(input$fev1)# Determina a cor com base na classificação
      
      div(
        style = paste0(
          "background-color: ", color, "; ",
          "color: white; ",
          "padding: 15px; ",
          "border-radius: 10px; ",
          "text-align: center; ",
          "font-weight: bold; ",
          "font-size: 18px;"
        ),
        current_gold_grade
      )
    })
    # Processa a classificação GOLD e cria uma visualização
    output$gold_classification <- renderText({ paste("Classificação GOLD:", current_gold_grade) })
    output$gold_classification_resumo <- renderText({ paste("Classificação GOLD:", current_gold_grade) })
    
    # Calcula a pontuação do CAT se o questionário for selecionado
    cat_score <- if (input$questionnaire_choice == "CAT") {
      sum(input$q1, input$q2, input$q3, input$q4, input$q5, input$q6, input$q7, input$q8)
    } else {
      NA
    }
    
    # Mostra resultados dinâmicos sobre eosinófilos e asma
    
    output$show_eosinophils_asthma <- reactive({
      total_exacerbations <- input$exacerbations_hosp + input$exacerbations_no_hosp
      total_exacerbations >= 2 || input$exacerbations_hosp >= 1
    })
    outputOptions(output, "show_eosinophils_asthma", suspendWhenHidden = FALSE)
    
    # Mostra resultados dinâmicos sobre dispneia (mMRC < 2 e poucas exacerbações)
    
    output$show_breathlessness <- reactive({
      total_exacerbations <- input$exacerbations_hosp + input$exacerbations_no_hosp
      total_exacerbations < 2 && input$mmrc < 2
    })
    outputOptions(output, "show_breathlessness", suspendWhenHidden = FALSE)
    
    # Trata inputs nulos ou NA de exacerbações de forma segura
    exacerbations_hosp <- ifelse(is.null(input$exacerbations_hosp) || is.na(input$exacerbations_hosp), NA, input$exacerbations_hosp)
    exacerbations_no_hosp <- ifelse(is.null(input$exacerbations_no_hosp) || is.na(input$exacerbations_no_hosp), NA, input$exacerbations_no_hosp)
    
    # Calcula o total de exacerbações
    total_exacerbations <- sum(c(exacerbations_hosp, exacerbations_no_hosp), na.rm = FALSE)
    
    # Define valores padrão para grupo e tratamento
    group <- "Por favor, insira todos os dados necessários"
    treatment <- "Por favor, insira todos os dados necessários"
    
    # Determina o grupo e o tratamento baseado nos critérios
    if ((!is.na(total_exacerbations) && total_exacerbations >= 2) || (!is.na(exacerbations_hosp) && exacerbations_hosp >= 1)) {
      # Grupo de alto risco
      group <- "Grupo E (Alto risco de exacerbação)"
      treatment <- if (!is.na(input$blood_eosinophils) && 
                       (input$blood_eosinophils >= 300 || input$asthma)) {
        "LABA + LAMA + ICS"  # Adiciona corticosteroide inalatório (ICS) em pacientes com eosinófilos elevados ou asma
      } else {
        "LABA + LAMA" # Combinação padrão de broncodilatadores
      }
    } else if (!is.null(input$mmrc) && input$mmrc >= 2 || 
               (!is.na(cat_score) && cat_score >= 10)) {
      # High symptom burden group
      group <- "Grupo B (Baixo risco, mas alto impacto de sintomas)"
      treatment <- "LABA + LAMA ou LABA ou LAMA"
    } else if (!is.null(input$mmrc) && input$mmrc < 2 && !is.na(total_exacerbations) && total_exacerbations < 2) {
      # Low risk and low symptom burden group
      group <- "Grupo A (Baixo risco e baixo impacto de sintomas)"
      treatment <- "Broncodilatador de longa duração"
    }
    
    # Renderiza os resultados do grupo e tratamento na interface
    output$group <- renderText({ paste("Grupo: ", group) })
    output$suggestion <- renderText({ paste("Tratamento Sugerido: ", treatment) })
    
    # Cria um resumo dinâmico com todas as informações do paciente
    output$summary <- renderText({
      resumo <- paste0(
        "Nome: ", input$name, "\n",
        "Idade: ", input$age, "\n",
        if (!is.na(cat_score)) paste("Pontuação CAT: ", cat_score, "\n"),
        if (input$questionnaire_choice == "mMRC") paste("Escala mMRC: ", input$mmrc, "\n"),
        "FEV1: ", input$fev1, "%\n",
        "Exacerbações com hospitalização: ", input$exacerbations_hosp, "\n",
        "Exacerbações sem hospitalização: ", input$exacerbations_no_hosp, "\n",
        if (group == "Grupo E") {
          paste("Eosinófilos: ", input$blood_eosinophils, " células/μL\n", 
                "Asma: ", ifelse(input$asthma, "Sim", "Não"), "\n")
        }
      )
      return(resumo)
    })
    
    # PDF Download Handler
    output$download_pdf <- downloadHandler(
      filename = function() {
        paste0("Relatorio_DPOC_", 
               gsub(" ", "_", input$name %||% "Paciente"), "_", 
               format(Sys.Date(), "%Y%m%d"), ".pdf")
      },
      content = function(file) {
        tryCatch({
          # Prepara os dados para o relatório
          report_data <- list(
            name = input$name %||% "Não Informado",
            age = input$age %||% "0",
            date = input$date %||% format(Sys.Date(), "%d de %B de %Y"),
            cc_number = input$cc_number %||% "Não Informado",
            utente_number = input$utente_number %||% "Não Informado",
            fev1 = input$fev1 %||% "0",
            gold_classification = gold_grade() %||% "Não Classificado",
            exacerbations_hosp = input$exacerbations_hosp %||% "0",
            exacerbations_no_hosp = input$exacerbations_no_hosp %||% "0",
            symptom_details = if (input$questionnaire_choice == "CAT") {
              sprintf("Pontuação CAT: %s", sum(input$q1, input$q2, input$q3, input$q4, input$q5, input$q6, input$q7, input$q8) %||% "0")
            } else {
              sprintf("Escala mMRC: %s", input$mmrc %||% "0")
            },
            treatment_suggestion = treatment %||% "Sem recomendação"
          )
          
          
          # Gerar PDF
          pdf_path <- generate_pdf_report(report_data)
          
          # Copiar PDF temporário para download
          if (file.exists(pdf_path)) {
            file.copy(pdf_path, file)
          } else {
            stop("Arquivo PDF não encontrado")
          }
        }, error = function(e) {
          # Mensagem de erro para o usuário
          showNotification(
            paste("Erro ao gerar PDF:", conditionMessage(e)), 
            type = "error"
          )
        })
       }
    )
  })
}
# Execute the app
shinyApp(ui, server)