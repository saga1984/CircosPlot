#!/usr/bin/env bash

#
# toma un transcriptoma ensamblado por trinity y genera archivos de entrada para circos plot
#

# menu de ayuda
AYUDA="
\n\tSINOPSIS: El script ejecuta modifica un archivo fasta de un trasncriptoma ensamblado por Trinity para circos plot
\n\tOPCIONES:
\n\t  OBLIGATORIOS
\t\t-i) ARCHIVO DE ENTRADA archivo fasta de nucleotidos (.fasta, .fa, .fna)
\n\t  OPCIONALES
\t\t-h) Muestra este MENU DE AYUDA
\n\tEJEMPLO: $(basename $0) -i Tinity.fa
"

# detener el script si no se introduce argumento (opcion)
if [[ $# -eq 0 ]]; then
   echo -e "${AYUDA}"
   exit 1
fi

modificar_archivo() {
   # generar archivo corto de gene_id y su longitud
   cat ${OPTARG} | grep ">" | cut -d " " -f "1,2" | sed "s/_i[0-9]//g" > gene_len_tmp.txt
   # crear inicio con chr[0-9] y eliminar >
   awk 'BEGIN {FS=OFS="_"} {split($3, a, /c/); print "chr" a[2] " " substr($0, 2)}' gene_len_tmp.txt > modified_gene_len_tmp.txt
   # ordenar por dege_id
   sort modified_gene_len_tmp.txt > sorted_modified_gene_len_tmp.txt
   # obtener conteos de inicio y final
   awk '{sum += substr($NF, 5); print $0, " " sum, " " prev_sum; prev_sum = sum }' sorted_modified_gene_len_tmp.txt > inicio_final_tmp.txt
   # eliminar la columna de len e invertir columnas inicio y final
   awk 'BEGIN {FS=" "; OFS="\t"} { print $1,$2,$5,$4}' inicio_final_tmp.txt > input_incomplete_tmp.tsv
   # agregar 0 en field 3 vacio en primera fila
   awk 'BEGIN {FS=OFS="\t"} {if ($3 ~ /^$/) $3 = 0; print $0 }' input_incomplete_tmp.tsv > preBED.tsv
}

tmp_files_remove() {
   # crear directorio de archivos temporales, solo si no existe
   if [[ ! -d tmp_files ]]; then
      mkdir tmp_files
   fi

   # mover archivos temporales a directorio correspondiente
   mv *_tmp.* tmp_files/

   # definir valor default de no eliminar archivos temporales
   DEFAULT="no"
   # mensaje preguntando si eliminar archivos temporales
   read -e -p "deseas remover archivos temporales (s/N):" REMOVE
   # comportamiento por defualt si se oprime boton ENTER (no borrar archivos temporales)
   REMOVE="${REMOVE:-${DEFAULT}}"
   # remover archivos temporales si se elige si
   if [[ "${REMOVE}" == "S" || "${REMOVE}" == "s" || "${REMOVE}" == "Si" || "${REMOVE}" == "si" ]]; then
      echo -e "\nRemoviendo carpeta de archivos temporales\n"
      rm -r tmp_files
      echo -e "\nModificacion de archivo terminada\n"
      exit 0
   # no remover archivos temporales si se elige no
   elif [[ "${REMOVE}" == "N" || "${REMOVE}" == "n" || "${REMOVE}" == "No" || "${REMOVE}" == "no" ]]; then
      echo -e "\nModificacion de archivo terminada\n"
      exit 0
   fi
}

# parsear argumentos
while getopts ":i:h" OPCIONES; do
   case ${OPCIONES} in
      i)
         # mensaje a stdout
         echo -e "\niniciando creacion de archivo guia para circos plot\n"
         # correr script que modifica fasta de transcriptomica para circos plot
         modificar_archivo
         # correr funcion para remover (o no) archivos temporales
         tmp_files_remove
         ;;
      h)
         # muestra menu de ayuda
         echo -e "${AYUDA}"
         exit 0
         ;;
      *)
         # parametro no valido
         echo -e "\n\tLa opcion -${OPTARG} no existe"
         echo -e "${AYUDA}"
         exit 1
         ;;
   esac
done
