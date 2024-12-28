import pandas as pd #usado para abrir o ficheiro
from tabulate import tabulate as tabu #usado para fazer as tabelas
import matplotlib.pyplot as plt #usado para fazer os gráficos
import seaborn as sns #usado para tornar os gráficos mais estéticos
from scipy import stats
from scipy.stats import chi2_contingency

ficheiro = "Sleep_health_and_lifestyle_dataset.csv"
df = pd.read_csv(ficheiro)  # Abrir/ler ficheiro

df['Sleep Disorder'] = df['Sleep Disorder'].fillna('None')
df['Sleep Disorder'] = pd.Categorical(
    df['Sleep Disorder'],
    categories=['None', 'Sleep Apnea', 'Insomnia'],
    ordered=True
)

bins = [64, 69, 74, 79, 84, 86]  # Limites para criar categorias/grupos
labels = ['65-69', '70-74', '75-79', '80-84',
          '85-86']  # Categorias para Heart Rate
# Criar a coluna com os grupos
df['Heart Rate'] = pd.cut(df['Heart Rate'], bins=bins,
                          labels=labels, right=True)

df['BMI Category'] = df['BMI Category'].replace('Normal', 'Normal Weight')

def categorizarPressao(leitura):
    sistolica = int(leitura.split('/')[0])

    if 110 <= sistolica < 120:
        return "110-119"
    elif 120 <= sistolica < 130:
        return "120-129"
    elif 130 <= sistolica < 140:
        return "130-139"
    else:
        return "140-149"

# Aplicando a função e criando a nova coluna 'BloodPressureGroups'
df['Blood Pressure Groups'] = df['Blood Pressure'].apply(categorizarPressao)

# resumo = df.describe(include='all')
# print(resumo)
#######################################################################################################


def analisar_variavel(df, variavel):
    # Criar figura com 3 subplots
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    serie = df[variavel]

    # Estatísticas e gráfico de frequência
    if pd.api.types.is_numeric_dtype(serie):  # Variável numérica
        estatisticas = {
            "Mean": serie.mean(),
            "Median": serie.median(),
            "Mode": serie.mode().iloc[0] if not serie.mode().empty else "N/A",
            "Standard Deviation": serie.std(),
            "Maximum": serie.max(),
            "Minimum": serie.min(),
            "Quantiles": serie.quantile([0.25, 0.5, 0.75]).to_dict(),
            "IQR": serie.quantile(0.75) - serie.quantile(0.25)
        }
        texto = "\n".join([f"{k}: {v}" for k, v in estatisticas.items()])
        axes[0].text(0.5, 0.5, texto, fontsize=11, ha="center", va="center")
        axes[0].axis("off")
        axes[0].set_title(f"Statistics: {variavel}", fontsize=14)
        sns.histplot(serie, kde=True, ax=axes[1], bins=20, color="skyblue")
        axes[1].set_title(f"Frequency Graphic: {variavel}", fontsize=14)

        # Boxplot para a relação com 'Sleep Disorder'
        sns.boxplot(hue='Sleep Disorder', y=variavel,
                    data=df, ax=axes[2], palette="pastel")
        axes[2].set_title(
            f"Distribution of {variavel} by Sleep Disorders", fontsize=14)

    else:  # Variável categórica

        freq_abs = serie.value_counts()
        freq_rel = serie.value_counts(normalize=True) * 100
        estatisticas = {
            "Mode": serie.mode().iloc[0] if not serie.mode().empty else "N/A",
            "Total Categories": serie.nunique(),
            "Absolute Frequencies": freq_abs.to_dict(),
            "Relative Frequencies(%)": freq_rel.round(2).to_dict(),
        }
        texto = "\n".join([f"{k}: {v}" if isinstance(
            v, (str, int, float)) else f"{k}:\n{v}" for k, v in estatisticas.items()])
        axes[0].text(0.5, 0.5, texto, fontsize=10,
                     ha="center", va="center", wrap=True)
        axes[0].axis("off")
        axes[0].set_title(f"Statistics: {variavel}", fontsize=14)

        sns.countplot(y=serie, ax=axes[1], palette="pastel", legend=False)
        axes[1].set_title(f"Frequencies Graphic: {variavel}", fontsize=14)

        # Gráfico de Barras Empilhadas
        if variavel != 'Sleep Disorder':
            barras = pd.crosstab(
                df[variavel], df['Sleep Disorder'], normalize='index')
            barras.plot(kind="bar", stacked=True,
                        ax=axes[2], colormap="Pastel1")
            axes[2].set_title(
                f"Distribution of {variavel} by Sleep Disorder", fontsize=14)
            axes[2].set_xlabel(variavel)
            axes[2].set_ylabel("Proportion")

    # Ajustar layout
    plt.tight_layout()
    plt.show()

# esta função usa a de cima. está terminada.


def mostrarVariaveis():
    listaVars = df.columns.tolist()
    while True:
        x = 1
        print("This menu allows the individual analysis of each variable and their association with the main variable: Sleep Disorders\nThe possible options for analysis are:")
        for i in range(1, len(listaVars)):
            print(x, " - ", listaVars[i])
            x = x+1
        # print(x, " - ", listaVars[len(listaVars)-1])
        # for i in range(1, len(listaVars)-2):
         #   print(x, " - ", listaVars[i])
          #  x = x+1
        # print(x, " - ", listaVars[len(listaVars)-1])
        print("99 - Go back to the main menu")
        resposta = (input("Write the number of the pretended option: "))
        print("")
        if resposta.isdigit():
            if int(resposta) > 0 and (int(resposta) < len(listaVars)):
                # if int(resposta) > 0 and (int(resposta) < len(listaVars)-2 or int(resposta) == 12):
                # if int(resposta) == 12:
                #    resposta = 13
                analisar_variavel(df, listaVars[int(resposta)])
            elif resposta == "99":
                break
            else:
                print(
                    "ERROR: Variable not found. Try again. \n-------------------||-------------------")
        else:
            print(
                "ERROR: Variable not found. Try again. \n-------------------||-------------------")


# esta função pretende mostrar correlações entre variáveis.
def mostraCorrs():
    while True:
        print("Welcome to the exploratory data analysis (EDA) between variables page.\nThe available options are: \n1 - All variables \n2 - Occupation vs Sleep disorder \n3 - Occupation vs Sleep Quality\n4 - Blood Pressure vs Sleep disorder\n5 - Stress level vs Sleep Disorder\n6 - Stress level vs Sleep Quality\n7 - Stress Level vs Sleep Duration\n99 - Go back to the main menu.")
        resposta = input("Your choice: ")
        if resposta == "1":
            # Matriz de correlações entre as variáveis numéricas
            df_numeric = df.select_dtypes(include=['number'])
            corr = df_numeric.corr().round(2)
            plt.figure(figsize=(25, 20))
            sns.heatmap(corr, annot=True, cmap='viridis')
            # Roda os nomes das variáveis para não ficarem sobrepostos
            plt.xticks(rotation=15, ha='right')
            plt.yticks(rotation=15, ha='right')

            plt.show()
        elif resposta == "2":
            # Occupation vs Sleep disorder
            sns.set_style("whitegrid")
            plt.figure(figsize=(10, 5))
            sns.histplot(data=df, y="Occupation", hue="Sleep Disorder",
                         multiple="stack", stat="percent", shrink=.8)

            plt.title(
                'Sample Distribution per Sleep Disorder and Occupation Graph 3')
            plt.xlabel('Percent')
            plt.show()
        elif resposta == "3":
            # Occupation vs Sleep Quality
            # Criar uma tabela de contingência
            contingency = pd.crosstab(
                df['Occupation'], df['Quality of Sleep'], normalize='index')

            # Gráfico de mosaico
            plt.figure(figsize=(10, 5))
            sns.heatmap(contingency, annot=True, cmap='YlGnBu', fmt='.2%')

            plt.title('Quality of Sleep per Occupation Graph 4')
            plt.xlabel('Quality of Sleep')
            plt.ylabel('Occupation')

            plt.tight_layout()
            plt.show()
        elif resposta == "4":
            # BP vs Sleep disorder
            sns.set_style("whitegrid")
            plt.figure(figsize=(10, 5))
            sns.histplot(data=df, y="Blood Pressure Groups", hue="Sleep Disorder",
                         multiple="stack", stat="percent", shrink=.8)

            plt.title(
                'Sample Distribution per Sleep Disorder and Sistolic Blood Pressure Graph 6')
            plt.ylabel('Sistolic Blood Pressure')
            plt.xlabel('Percent')
            plt.show()
        elif resposta == "5":
            # stress level vs sleep disorder
            # Criar uma tabela de contingência
            contingency = pd.crosstab(
                df['Stress Level'], df['Sleep Disorder'], normalize='index')

            # Gráfico de mosaico
            plt.figure(figsize=(10, 5))
            sns.heatmap(contingency, annot=True, cmap='YlGnBu', fmt='.2%')
            plt.title('Sleep Disorders per Stress Levels Graph 9')
            plt.xlabel('Sleep Disorders')
            plt.ylabel('Stress Level')
            plt.tight_layout()
            plt.show()
        elif resposta == "6":
            # Stress level vs Sleep Quality
            # Criar uma tabela de contingência
            contingency = pd.crosstab(
                df['Stress Level'], df['Quality of Sleep'], normalize='index')

            # Gráfico de mosaico
            plt.figure(figsize=(10, 5))
            sns.heatmap(contingency, annot=True, cmap='YlGnBu', fmt='.2%')
            plt.title('Quality of Sleep per Stress Level Graph 10')
            plt.xlabel('Quality of Sleep')
            plt.ylabel('Stress Level')
            plt.tight_layout()
            plt.show()
        elif resposta == "7":
            # Stress Level e Sleep Duration
            plt.figure(figsize=(10, 5))
            sns.boxplot(x='Stress Level', y='Sleep Duration',
                        palette='deep', data=df)
            plt.title(
                'Sample Distribution of Sleep Duration per Stress Level Graph 11')
            plt.xlabel('Stress Level')
            plt.ylabel('Sleep Duration (hours)')
            plt.show()
        elif resposta == "99":
            break
        else:
            print("Unknown answer. Please, try again.\n")

def QuiSquareTest_bmi():
    from scipy.stats import chi2_contingency
    # Criar a tabela cruzada
    contingency_table = pd.crosstab(df['BMI Category'], df['Sleep Disorder'])

    chi2, p, dof, expected = chi2_contingency(contingency_table)
    print(f"Qui-Square Test between BMI Categories and Sleep Disorders:")
    print(f"Qui-Square Statistics: {chi2}")
    print(f"p-value: {p}")
    if p < 0.05:
        print("There is a significant association between BMI Categories and Sleep Disorders.")
    else:
        print("There is not a significant association between BMI Categories and Sleep Disorders.")

def QuiSquareTest_sleep():
    cat = pd.cut(df['Quality of Sleep'], bins=[0, 4, 7, 10], labels=['Low', 'Medium', 'High'])
    contingency_table2 = pd.crosstab(cat, df['Sleep Disorder'])
    chi2_result = stats.chi2_contingency(contingency_table2)
    print(f"Chi-Square Statistic: {chi2_result.statistic:.2f}")
    print(f"P-value: {chi2_result.pvalue:.2e}")
    if chi2_result.pvalue < 0.05:
        print("There is a significant association between sleep quality and sleep disorders.")
    else:
        print("No significant association between sleep quality and sleep disorders.")

def datasetStr():
    # ESTATISTICA DESCRITIVA
    print("Showing the first 5 observations...")
    print(df.head())
    # print(df.columns)

    shape = df.shape
    print('The dataset has {} observations and {} variables'.format(
        shape[0], shape[1]))

    # to check NAs and variable type
    # df.info()
    input("Press enter to continue and see the continuos variables.")
    # vars numericas
    estatistica1 = df.describe().round(2).drop('count')
    print(tabu(estatistica1, headers='keys', tablefmt='fancy_grid'))
    input("Press enter to continue and see the categorical variables.")

    # vars categoricas
    estatistica2 = df.describe(include='O').drop(
        'count').drop('Blood Pressure', axis=1)
    estatistica2['Sleep Disorder'] = df['Sleep Disorder'].describe()
    print(tabu(estatistica2, headers='keys', tablefmt='fancy_grid'))

    # percentagem de doentes
    print("Sleep disorders proportions in this sample:")
    percent = df['Sleep Disorder'].value_counts(normalize=True) * 100
    percent = percent.round(2)
    print(percent)


def saveOnFile():
    BMI_order = df['BMI Category'].value_counts().index

    # Abrir o arquivo em modo de escrita (irá sobrescrever ou criar um novo arquivo)
    with open('Dados_e_gráficos_Projeto.txt', 'w', encoding="utf-8") as file:
        # Escrever título e primeiras observações (head)
        file.write('### First 5 Observations of the DataFrame ###\n')
        file.write(str(df.head()))  # Salva as primeiras 5 linhas
        file.write('\n\n')  # Espaço entre seções

    # Escrever as variáveis do DataFrame
        file.write('### Variables of the DataFrame ###\n')
        file.write(str(df.columns))  # Salva os nomes das colunas
        file.write('\n\n')  # Espaço entre seções

    # Teste Qui-Quadrado (se aplicável)
        from scipy.stats import chi2_contingency

        contingency_table = pd.crosstab(
            df['BMI Category'], df['Sleep Disorder'])
        chi2, p, dof, expected = chi2_contingency(contingency_table)

        file.write('### Qui-Square Test between BMI Categories and Sleep Disorders ###\n')
        file.write(f"Qui-Square Statistics: {chi2}\n")
        file.write(f"p-value: {p}\n")
        if p < 0.05:
            file.write(
                "There is a significant association between BMI Categories and Sleep Disorders.\n")
        else:
            file.write(
                "There is not a significant association between BMI Categories and Sleep Disorders.\n")
        file.write('\n\n')

        # Distribuição das observações pela BMI
        plt.subplot(1, 2, 1)
        plt.gca().set_title('BMI Category Distribution', fontweight='bold')
        sns.countplot(x='BMI Category', hue='BMI Category',
                      palette='Set2', data=df, order=BMI_order)
        plt.xlabel('BMI Category', fontweight='bold', fontsize=12)
        plt.ylabel('Frequency', fontweight='bold', fontsize=12)
        plt.xticks(rotation=0)
        # Distribuição das observações pela sleep disorder
        plt.subplot(1, 2, 2)
        plt.gca().set_title('Sleep Disorder Distribution', fontweight='bold')
        sns.countplot(x='Sleep Disorder', hue='Sleep Disorder',
                      palette='Set2', data=df)
        plt.xlabel('Sleep Disorder', fontweight='bold', fontsize=12)
        plt.ylabel('Frequency', fontweight='bold', fontsize=12)
        # Salvar o gráfico como imagem
        plt.savefig('BMI_Distribution_and_SleepDisorder_Distribution.jpeg')
        plt.close()  # Fechar o gráfico
        file.write(
            "### Graphic: Distribution of BMI Categories and Sleep Disorders ###\n")
        file.write(
            "Graphic saved as 'BMI_Distribution_and_SleepDisorder_Distribution.jpeg'.\n\n")

        # Gráfico BMI vs Sleep Duration
        plt.figure(figsize=(8, 6))
        BMI_order = df['BMI Category'].value_counts().index
        sns.boxplot(x='BMI Category', y='Sleep Duration',
                    data=df, palette='Set2', order=BMI_order)
        plt.title('Sleep Duration vs BMI Category', fontweight='bold')
        plt.xlabel('BMI Category', fontweight='bold', fontsize=12)
        plt.ylabel('Sleep Duration (hours)', fontweight='bold', fontsize=12)
        plt.savefig('BMI_vs_Sleep_Duration.jpeg')
        plt.close()
        file.write("### Graphic: BMI Categories vs Sleep Duration ###\n")
        file.write("Graphic saved as 'BMI_vs_Sleep_Duration.jpeg'.\n\n")
        # Gráfico BMI vs Quality of Sleep
        plt.figure(figsize=(8, 6))
        sns.boxplot(x='BMI Category', y='Quality of Sleep',
                    data=df, palette='Set2', order=BMI_order)
        plt.title('Sleep Quality Distribution by BMI Category',
                  fontweight='bold')
        plt.xlabel('BMI Category', fontweight='bold', fontsize=12)
        plt.ylabel('Quality of Sleep (0 to 10)',
                   fontweight='bold', fontsize=12)
        plt.savefig('BMI_vs_Sleep_Quality.jpeg')
        plt.close()
        file.write("### Graphic: BMI Categories vs Quality of Sleep ###\n")
        file.write("Graphic saved as 'BMI_vs_Sleep_Quality.jpeg'.\n\n")


        # Sleep duration vs Sleep disorder
        plt.figure(figsize=(12, 6))
        plt.subplot(1, 2, 1)
        sns.boxplot(x='Sleep Disorder', y='Sleep Duration',
                    data=df, palette='coolwarm')
        plt.title('Sleep Duration vs Sleep Disorder', fontweight='bold')
        plt.xlabel('Sleep Disorder', fontweight='bold', fontsize=12)
        plt.ylabel('Sleep Duration (hours)', fontweight='bold', fontsize=12)

# esta é a main(). detem o UI


def main():
    print("\n\nExploring Sleep Health and Lifestyle Patterns Through Python Data Analysis \nWelcome to our project developed at the UC of Programming and Algorithms in Sciences. \nThe project was made by \nAfonso Carreira, nº mec. 107988 \nBeatriz Lavado, nº mec. 102495 \nTomás Geraldes, nº mec. 107508 \n-------------------||------------------- \n\nThe dataset used for this work was 'Sleep_health_and_lifestyle_dataset.csv' \nTo get this dataset, visit https://www.kaggle.com/datasets/uom190346a/sleep-health-and-lifestyle-dataset/data \n")
    while True:
        print("-------------------||-------------------")
        print("Sleep is a physiological need present in animals, but it is increasingly being neglected, particularly due to work pressures, heavy workloads, and productivity demands. In summary, this study aims to show how lifestyles can be associated with sleep quality and vice versa, using Python code to visualize the results and conclusions of the work.")
        print(
            "To learn more about the work, the following options are available: ")
        # fazer contxtualização sobre o sono
        print("1 - Contextualization about sleep.")
        # Ver os nomes, ver medias, medianas, sd, etc; ver gráficos de frequencias
        print("2 - Study variables. (see graphics and dispersion measures)")
        print("3 - Most relevant variables' associations") 
        print("4 - Qui-Square Test between BMI Categories and Sleep Disorders") 
        print("5 - Qui-Square Test between Sleep Quality and Sleep Disorders")
        print("6 - See dataset structure")
        print("7 - Save data on a file")
        print("8 - See credits")
        print("99 - Terminate")
        resposta = input("Your answer: ")
        print("\n-------------------||-------------------")
        if resposta == "1":
            print("Sleep is an essential biological process that allows the body and mind to recover energy and function properly. \nDuring sleep, important processes occur, such as memory consolidation, tissue repair, immune system strengthening, and hormone regulation. \nSleep is divided into cycles with different phases, like REM and non-REM sleep, each playing a specific role in physical and mental health. \nSleep quality and duration are influenced by internal factors, such as the circadian rhythm, and external ones, like lifestyle, stress levels, and health conditions. \nPoor or insufficient sleep is linked to issues like fatigue, difficulty concentrating, mood changes, and an increased risk of chronic diseases, highlighting its importance for overall well-being.\n")
        elif resposta == "2":
            mostrarVariaveis()
        elif resposta == "3":
            mostraCorrs()
        elif resposta == "4":
            QuiSquareTest_bmi()
        elif resposta == "5":
            QuiSquareTest_sleep()
        elif resposta == "6":
            datasetStr()
        elif resposta == "7":
            saveOnFile()
        elif resposta == "8":
            print("Project developed at the UC of Programming and Algorithms in Sciences \nThe students:\nAfonso Carreira, 107988 \nBeatriz Lavado, 102495 \nTomás Geraldes, 107508 \n \n Obrigado aos nossos pais por suportarem a nossa estadia na Universidade.")
        elif resposta == "99":
            break
        else:
            print("ERROR: Unknown answer. Try again. \n")


main()