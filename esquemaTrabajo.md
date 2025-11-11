# 1. Seleccionas una ruta biológica
Por ejemplo, “Apoptosis humana” de KEGG o Reactome. → obtienes su lista de genes (~100 o menos).

# 2. Tomas esa lista como conjunto de semillas (seed genes). 
La idea es esos genes son los “miembros conocidos”.

# 3. Propagas esa información sobre una red PPI (p. ej. STRING, BioGRID, o una versión recortada).
- Usamos DIAMOND o GUILD (es posible usar Random Walk, pero mejor usar estos que ya hemos usado)
- Obtienes una puntuación de “proximidad” para todos los genes de la red. 

# 4. Ordenas por puntuación y te quedas con los nuevos genes top-k que no estaban en la ruta original. 
Estos son tus candidatos a formar parte de la ruta. (en el caso de diamond, filtramos por menor p-value o rank, en guild si se obtiene una puntuacion continua de propagacion)

PASO ADICIONAL: analizar diferentes parametros de la red con la libreria networkx, como modularidad y otros

# 5. Análisis funcional
- Enriquecimiento GOKEGG de los genes originales → “control”.
- Enriquecimiento de los nuevos genes → comparas si salen los mismos términos (“coherencia”) o términos complementarios (“extensión funcional”).
- A veces, los nuevos genes completan subprocesos (“regulación positiva”, “mantenimiento mitocondrial”, etc.).

# 6. Conclusión biológica
- “Estos genes vecinos podrían ser participantes adicionales en la vía de apoptosis.”
- “La propagación resalta genes que interaccionan con las caspasas, coherentes con el proceso.”

# 7. Visualización
- Subred en Cytoscape o Python (pyvisnetworkx)
	- Azul genes originales
	- Rojo genes añadidos tras propagación
	- Tamaño puntuación o grado


**NO OLVIDAR HACER LOS .SH PARA QUE SE EJECUTEN TODOS LOS SCRIPTS NECESARIOS**  