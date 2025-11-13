# 1. Selección de la ruta biológica
Para este proyecto, hemos seleccionado la ruta **"Autophagy - animal" (KEGG hsa04140)**. A partir de la base de datos KEGG o Reactome se obtiene la lista de genes que participan en la vía canónica, que es la parte funcional del proceso. 

# 2. Conjunto de semillas (seed genes)
Estos genes de la vía autofágica serán los que se van a usar como **semillas** en la red de interacción proteína-proteína. La idea es que las semillas representen los miembros conocidos y validados experimentalmente del proceso. A partir de estos miembros se intentarán identificar **nuevos genes potencialmente relacionados** mediante la propagación de información a través de la red.

# 3. Propagación de la información sobre la red PPI
Para la propagación, se utilizará una red de alta cobertura, como **STRING** o **BioGRID** y se aplica un algoritmo de propagación, como **DIAMOnD** o **GUILD**, que permiten priorizar genes en función de su proximidad topológica a las semillas.

Cada gen obtendrá una puntuación de cercanía respecto al conjunto inicial de genes de autofagia. De esta forma podremos ordenar todos los genes de la red según su relación potencial con el proceso autofágico.

# 4. Genes candidatos
Una vez obtenida la red según la puntuación, seleccionamos los genes que no formen parte de la vía original de autofagia y que tengan una posición alta en el ranking (p.ej. los 100 primeros).

Estos genes serán los candidatos a estar implicados en la autofagia, aunque no figuren en las bases de datos de rutas biológicas. Posteriormente se pueden dividir por rangos de puntuación, para analizar cómo varía la coherencia funcional con la distancia.

# 5. Análisis estructural y funcional
Primero se deben estudiar las propiedades estructurales de la subred formada por los genes de autofagia y los nuevos candidatos. Podemos calcular distintos parámetros para ello, como el **grado**. o la **modularidad**. Esto nos va a ayudar a comprender cómo se integran los nuevos genes dentro del contexto de la autofagia.

Para el enriquecimiento funcional, debemos separarlo en: 

* Enriquecimiento GOKEGG de los genes originales: “Control”.
* Enriquecimiento de los nuevos genes: Comparar si salen los mismos términos (“coherencia”) o términos complementarios (“extensión funcional”).

A veces, los nuevos genes completan subprocesos (“regulación positiva”, “mantenimiento mitocondrial”, etc.).

# 6. Conclusión biológica
Ejemplos de conclusiones a las que podríamos llegar:

* *Estos genes vecinos podrían ser participantes adicionales en la vía de autofagia.*
* *La propagación resalta genes que interaccionan con las caspasas, coherentes con el proceso.*

# 7. Visualización
- Subred en Cytoscape o Python (pyvisnetworkx)
	- Azul genes originales
	- Rojo genes añadidos tras propagación
	- Tamaño puntuación o grado


**NO OLVIDAR HACER LOS .SH PARA QUE SE EJECUTEN TODOS LOS SCRIPTS NECESARIOS**  