# Proyecto de Análisis Funcional y Propagación en Redes

## 1. Introducción

Durante este proyecto se ha desarrollado un pipeline integrado para realizar la descarga de los genes de la autofagia, la **propagación de la red con DIAMOND** y el **análisis computacional** de genes relacionados, combinando dos enfoques:

* **Enriquecimiento funcional (GSEA/Enrichr)**. Esto lo utilizamos para identificar rutas biológicas, procesos celulares y mecanismos moelculares asociados tanto a los **genes semilla** como a los **genes candidatos**.
* **Análisis estructural/topológico de redes PPI**. Este análisis se realiza para estudiar el comportamiento de estos genes dentro de una **red de interacciones proteicas**, identificando nodos clave, genes puente, modularidad y roles dentro del propio sistema celular.

El objetivo de este proyecto ha sido expandir la lista de genes de una ruta biológica de interés, que en este caso ha sido **"Autophagy - animal" (KEGG hsa04140)**. Así pues, se va a validar **si los genes candidatos podrían estar funcionalmente o estructuralmente relacionados con la autofagia** y priorizar aquellos con mayor relevancia.

## 2. Estructura del repositorio

HAB_Proyecto/
├── launch.sh                   # Lanzador principal del pipeline secuencial
├── requirements.txt            # Dependencias de Python.
├── scripts/
│   ├── descargar_ruta.py       # Descarga de genes semilla (PASO 1)
│   ├── propagacion_diamond.py  # Ejecuta el algoritmo DIAMOnD (PASO 2).
│   ├── enriquecimiento_funcional.py             # Análisis funcional y estructural de resultados (PASO 3).
│   └── descarga_red_string.py  # Script adicicional: descarga la red de STRINGDB
├── data/                     
│   ├── genes_ruta.txt          # Lista de genes semilla (Output PASO 1).
│   └── string_network_filtered_hugo-400.tsv # Red PPI filtrada y mapeada a HUGO.
└── results/     
    ├── analisis_estructural.tsv     
    ├── comparacion_enriquecimiento.txt 
    ├── connected_seed_genes.tsv     
    ├── isolated_seed_genes.tsv     
    ├── diamond_added_genes.tsv    
    ├── diamond_results.tsv        
    ├── enriquecimiento_candidatos.tsv 
    ├── enriquecimiento_semillas.tsv 
    ├── subred_autofagia.graphml    
    ├── subred_enriquecida.graphml   
    ├── diamond_network.png          
    └── .png                    # Otros gráficos generados del análisis funcional               

## 3. Instrucciones de ejecución

Primero de todo, hay que instalar las dependencias necesarias para python. Asegúrate de estar en el directorio raíz del proyecto (HAB_Proyecto/) y ejecuta:
```bash
pip install -r requirements.txt
```

Después, el flujo de trabajo completo se lanza con un script de shell llamado *launch.sh*, utilizando un terminal compatible con Unix (Linux, macOS, o Git Bash/WSL en Windows) con el siguiente comando:
```bash
./launch.sh
```

## 4. Metodología y análisis
### 4.1. Descargar los genes la ruta

El script **descargar_ruta.py** se encarga de descargar automáticamente los genes de la autofagia desde las librerías de Enrichr y guardarlos en un archivo de texto.

En este proyecto lo usaremos para:

* Obtener los genes semilla de la autofagia directamente de una base de datos funcional.
* Generar el archivo `data/genes_ruta.txt`, que luego se usará como entrada en el script de **propagacion_diamond.py**

#### Funcionamiento del script

1. Selecciona una librería funcional de Enrichr.
2. Busca una ruta cuyo nombre contenga una palabra clave.
3. Descarga todos los genes asociados a esa ruta.
4. Localiza la primera ruta cuyo nombre contiene "autophagy", y guarda:
    * El nombre completo de la ruta encontrada.
    * La lista de genes asociados.
5. Escribe los genes en un archivo de texto `data/genes_ruta.txt`, un gen por línea en formato HUGO.

#### Resultados generados
**1. data/genes_ruta.txt**

Contiene la lista de genes de la ruta de la autofagia, un gen por línea.

### 4.2. Propagación DIAMOnD

En el sript de **propagacion_diamond.py** implementamos el algoritmo **DIAMOnD** sobre una red PPI ya proporcionada, en este caso de STRINGDB filtrada con score superior a 400. Se ha seleccionado este algortimo de propagación ya que se trata de una herramienta con gran robustez estadística que se enfoca específicamente en la detección de módulos de genes, es decir, busca el camino que maximiza la coherencia interna del grupo, lo que resulta en un módulo funcional más cohesivo y biológicamente interpretable. Este script tiene el propósito de:

* Tomar como entrada un conjunto de **genes semilla** (HUGO symbols) relacionados con autofagia.
* Propagar sobre una red de interacción (STRING) que se proporciona ya de inicio en la carpeta **data** y se añaden a nuestros genes semillas, otros nodos (genes) que estén significativamente conectados con el cluster de semillas. Para garantizar robustez, hemos filtrado todos aquellos genes con un umbral superior a 700.
* Separar claramente:
    * Genes semilla **conectados** a la red filtrada.
    * Genes semilla **aislados**.
    * Genes **añadidos por DIAMOnD** como candidatos.

Es el primer paso del pipeline: genera los genes candidatos y define qué semillas se pueden se pueden usar realmente en la red.

#### Funcionamiento del script

1. Carga de genes semilla desde un archivo de texto.
2. Carga de la red PPI y filtrado por score.
3. Construcción del grafo con NetworkX.
4. Comprobación de qué genes semilla están presentes y conectados en la red.
5. Ejecución de **DIAMOnD**:
    * A partir de las semillas válidas.
    * Añadiendo un número definido de nodos mientras existen candidatos.
6. Guardado de resultados:
    * Genes semilla conectados.
    * Genes semilla aislados.
    * Genes candidatos.
7. Visualización de la red enriquecida en formato imagen.

#### Resultados generados

**1. diamond_results.tsv**

Contiene lso genes añadidos por DIAMOnD, es decir, los genes candidatos.

**2. connected_seed_genes.tsv**

Genes semilla que están presentes y conectados en la red PPI filtrada.

**3. isolated_seed_genes.tsv**

Genes semilla que no están conectados en la red a partir del umbral score. Estos genes no participan en la propagación y pueden comentarse como limitación de la red.

**4. diamond_network.png**

Grafo de la subred enriquecida que incluye; nodos azules (genes semilla conectados) y nodos naranjas (genes añadidos por DIAMOnD). Permite una visualización rápida de cómo se expande el cluster de autofagia en la red PPI.

#### Interpretación de resultados

![Network](/results/diamond_network.png)

La mayoría de los genes semilla aparecen **conectados dentro de la red STRING**, formando un clúster denso en el centro del grafo. Esto nos indica que la autofagia es un módulo altamente cohesivo, con muchas interacciones internas bien documentadas.

Los genes naranjas corresponden a los candidatos generados por el algoritmo. Su patron general es:

* Muchos candidatos actúan como **puentes** entre los genes de autofagia y otros módulos celulares.
* Otros aparecen en ramas periféricas, lo que sugiere roles **indirectos o reguladores**.

Hay dos genes que aparecen como no conectados en la red filtrada, que al no tener suficientes interacciones con score > 700 **no entran en la propagacion DIAMOnD** y quedan fuera del análisis. Podrían ser genes que actúan de manera muy periférica o específica en la ruta de autofagia y por ello no encontramos interacciones.

### 4.3. Enriquecimiento funcional

El script **enriquecimiento_funcional.py** realiza un análisis sobre dos grupos de genes:

* **Genes semilla** -> Son genes ya conocidos de la ruta de la autofagia.
* **Genes candidatos** -> Son genes derivados del análisis predictivo.

Para ambos conjuntos se ejecuta un enriquecimiento con:

* **KEGG 2021 hUMAN**
* **GO Biological Process 2021**

Después se comparan los resultados entre sí y se generan gráficos interpretables.

#### Funcionamiento del script

1. Carga listas de genes desde archivos TSV.
2. Realiza enriquecimiento funcional para cada conjunto (semillas y candidatos).
3. Compara los 20 términos más significcativos entre ambos grupos.
4. Genera visualizaciones (-log10(FDR)).
5. Guarda tablas completas de resultados.
6. Guarda un informe comparativo con términos comunes y exclusivos.

#### Resultados generados

**1. enriquecimiento_semillas.tsv**

Tabla completa con rutas de autofagia y procesos relacionados.

**2. enriquecimiento_candidatos.tsv**

Tabla con rutas enriquecidas en los genes predichos.

**3. comparacion_enriquecimiento.txt**

Informe con términos comunes, términos exclusiovos de semillas y términos exclusivos de candidatos.

**4. Gráficos en PNG**

* `top_terms_Semillas.png`
* `top_terms_Candidatos.png`

Permiten identificar visualmente los términos más significativos (-log10(FDR)).

#### Interpretación de resultados

Los resultados muestran un patrón claro y biológicamente coherente entre lso dos grupos analizados.

##### Enriquecimiento de genes semilla

![Semillas](results\top_terms_Semillas.png)

Los genes semilla recuperan exactamente las rutas esperadas del proceso autofagico. La presencia simultánea de reguladores como **mTOR**, **AMPK** y **FoxO** refuerza que el conjunto está bien definido y sirve como un control positivo del análisis.

Esto sirve de validación de que el pipeline de enriquecimiento funciona correctamente.

##### Enriquecimiento de genes candidatos

![Candidatos](results\top_terms_Candidatos.png)

Las rutas que se muestran no son directamente de la autofagia, sino que son procesos de **señalización upstream** que regulan o modulan la autofagia.

En concreto, Rap1, Ras y MAPK son módulos clave que controla la señalización del estrés, la prloferación celular, la regulación metabólica y la activación de vías como mTOR o AMPK.

Esto nos sugiere que lso genes candidatos podrían actuar antes que lso genes semilla, podrían aprticipar en la activación o represión del proceso autofágico o podrían integrar señales externas en la maquinaria autofágica.

La presencia de procesos de neurodegeneración también es consistente, ya que la autofagia está fuertemente implicada en la degradación de proteínas agregadas en estas patologías.

#### Valor biológico del análisis

El patrón observado es el esperado en un análisis de expansión de redes:

* Las semillas confirman el proceso de interés.
* Los candidatos revelan **reguladores**, **moduladores** y **contextos patológicos** relacionados.
* El solapamiento pequeño pero coherente sugiere que los candidatos nos están directamente en el núcleo de la vía, sino arriba, conectando la autofagia con señalización celular y enfermedades asociadas.

## 5. Conclusión final 

El pipeline completo implementado nos proporciona una visión integradora sobre el papel potencial de nuevos genes en la dinámica de la autofagia.

En primer lugar, los **genes semilla** recuperados desde KEGG conforman un conjunto biológicamente robusto y coherente, lo que se confirma por su conectividad interna y por el enriquecimiento en procesos de autofagia, señalización lisosomal y regulación metabólica. Esto estableme un marco muy fiables sobre el que construir la propagación de la red.

La **propagación mediante DIAMOnD** añade un conjunto de **genes candidatos** que nos muestran patrones topológicos consistentes con roles biológicos relevantes. Muchos de ellos son **ndoods puente** entre el módulo de la autofagio y otros sistemas celulares, mientras que otros expanden el clúster hacia otras regiones periféricas. La exclusión de algunos genes semilla por falta de conectividad también muestra la importancia del filtrado por score y las limitaciones inherentes a las bases de datos PPI.

El **enriquecimiento funcional** refuerza esta interpretación:

* Los **genes semilla** reconstituyen fielmente la ruta de autofagia y sus reguladores inmediatos.
* Los **genes candidatos** se asocian a rutas upstream, así como a procesos de neurodegeneración y señalización de estrés, todos ellos moduladores conocidos del inicio y la intensidad del proceso autofágico.

Por su parte, el **análisis estructural/topológico** confirma que algunos de estos candidatos ocupan posiciones estratégicas en la red PPI, destacando por su grado, centralidad o betweenness, propiedades que los señalan como hubs reguladores o conectores clave entre módulos funcionales.

En conjunto, los resultados apoyan que los genes identificados por DIAMOnD no son simples vecinos en la red, sino **posibles reguladores de la autofagia**, actuando tanto a nivel de señalización como de integración de respuestas celulares. Este pipeline permite priorizar de forma objetiva un subconjunto de genes candidatos con mayor relevancia biológica y justificas su estudio posterior, ya sea mediante análisis computacionales adicionales o validación experimental en laboratorio.

## 6. Anexo

Se ha incluido en nuestra carpeta de scripts, el archivo **descarga_ruta_string.py** para mostrar el proceso que se tendría que realizar para no usar una red ya proporcionada. Se ha decidido omitir este paso ya que la descarga puede resultar lenta y puede dar resultados no reproducibles si se actualiza la base de datos de **STRINGDB**.

## 7. Bibliografía

Ghiassian, S. D., Menche, J., & Barabási, A. L. (2015). A disease module is a set of proteins with altered connectivity in disease. Nature Communications, 6(1), 1-13.

Levine, B., & Kroemer, G. (2019). Autophagy in the pathogenesis of disease. Cell, 176(1-2), 24-42.

Kuleshov, M. V., et al. (2016). Enrichr: a comprehensive gene set enrichment analysis resource. Nucleic Acids Research, 44(W1), W90-W97.