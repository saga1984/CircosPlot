#' Circos plot de DEGs usando coordenadas tipo preBED (Trinity)
#'
#' Lee una tabla de DEGs y un TSV con coordenadas estilo preBED (derivado de un
#' FASTA ensamblado con Trinity) y genera dos *circos plots* (Up y Down).
#' Exporta además archivos CSV con las posiciones y valores de \code{log2FC}.
#'
#' @param ruta_degs Directorio donde se encuentra el CSV de DEGs.
#' @param archivo_degs Nombre del archivo CSV de DEGs. Debe contener al menos
#'   las columnas \code{Gene_id}, \code{DEGs} (valores "Upregulated"/"Downregulated")
#'   y \code{log2FC}. El prefijo antes del primer guion bajo (\_) se usa como
#'   nombre de tratamiento para el subdirectorio de salida.
#' @param ruta_transcriptoma_modificado Directorio del TSV preBED.
#' @param archivo_transcriptoma_modificado Nombre del TSV preBED con 4 columnas:
#'   \code{Chr}, \code{Gene_id}, \code{Start}, \code{End} (separado por tabulador).
#' @param deg_col_up Color para los genes \emph{Upregulated} (puntos/líneas).
#' @param deg_col_down Color para los genes \emph{Downregulated} (puntos/líneas).
#' @param formatos Vector de formatos/gráficos a exportar (p. ej. \code{"tiff"},
#'   \code{"jpeg"}). Cada formato abre un dispositivo gráfico con \code{match.fun}.
#' @param resolucion Resolución en ppp (solo aplica a formatos ráster como TIFF/JPEG/PNG).
#'
#' @details
#' Se crea un subdirectorio \code{<TRATAMIENTO>_circosplot} dentro de
#' \code{ruta_transcriptoma_modificado}. Allí se guardan:
#' \itemize{
#'   \item \code{circos_up_source.csv} y \code{circos_down_source.csv}:
#'         tablas completas tras el \code{merge} con DEGs.
#'   \item \code{circos_up.csv} y \code{circos_down.csv}:
#'         columnas limpias \code{Chr, Inicio, Fin, log2FC}.
#'   \item \code{circos_up.<formato>} y \code{circos_down.<formato>}: figuras.
#' }
#' Los cromosomas se colorean con \code{circlize::rand_color(1)}. Para
#' reproducibilidad, establece una semilla antes de llamar a la función.
#'
#' @return \code{NULL} (efectos secundarios: crea directorios/archivos y guarda figuras).
#'
#' @examples
#' \dontrun{
#' circosplot_degs(
#'   ruta_degs = "2024/posdoc_CIBA/degs/",
#'   archivo_degs = "NACL_degs.csv",
#'   ruta_transcriptoma_modificado = "2024/posdoc_CIBA/",
#'   archivo_transcriptoma_modificado = "preBED.tsv",
#'   formatos = c("tiff","jpeg"),
#'   deg_col_up = "#0073C2FF",
#'   deg_col_down = "#bb0c00",
#'   resolucion = 600
#' )
#' }
#'
#' @seealso \href{https://cran.r-project.org/package=circlize}{circlize}
#' @keywords circos transcriptomica DEGs visualizacion
#' @import circlize
#' @importFrom dplyr group_by summarise
#' @importFrom utils read.csv write.csv
#' @importFrom grDevices dev.off
#' @author Valentín
#' @export


# funcion
circosplot_degs <- function(
    ruta_degs,
    archivo_degs,
    ruta_transcriptoma_modificado,
    archivo_transcriptoma_modificado,
    deg_col_up = "blue",
    deg_col_down = "red",
    formatos = "jpeg",
    resolucion = 300) {
  
  # programa a instalar
  programa = "circlize"
  
  # 
  if(programa %in% installed.packages() == FALSE){
    install.packages(programa)
  } else {
    print(paste("el programa", programa, "ya esta instalado", sep = " "))
  }
  
  library(circlize) 
  library(tidyverse)
  
  # tratamiento
  tratamiento <- gsub("\\_.*", "", archivo_degs)
  
  # importar archivos
  global <- read.csv(paste(ruta_transcriptoma_modificado, 
                           archivo_transcriptoma_modificado, 
                           sep = ""),
                     strip.white = TRUE,
                     header = FALSE,
                     sep = "\t")
  
  # definir nombre de columnas
  colnames(global) <- c("Chr", "Gene_id", "Start", "End")
  
  # crear data frames de degs con nombres cortos
  degs <- read.csv(paste(ruta_degs, archivo_degs, sep = ""), 
                   sep = ",",
                   strip.white = TRUE)
  
  # filtrar para remover not significant genes up y down
  degs_up <- degs[degs$DEGs == "Upregulated",]
  degs_down <- degs[degs$DEGs == "Downregulated",]
  
  # filtrar global en base a degs up y down
  global_up <- global[global$Gene_id %in% degs_up$Gene_id,]
  global_down <- global[global$Gene_id %in% degs_down$Gene_id,]
  
  # obtener valores promedio por Gene_id inducidos
  global_up_group <- global_up %>% 
    group_by(Chr, Gene_id) %>%
    dplyr::summarise(Inicio = round(mean(Start)), 
                     Fin = round(mean(End, na.rm = TRUE)))
  
  # obtener valores promedio por Gene_id reprimidos
  global_down_group <- global_down %>% 
    group_by(Chr, Gene_id) %>%
    dplyr::summarise(Inicio = round(mean(Start)), 
                     Fin = round(mean(End, na.rm = TRUE)))
  
  # agregar columna de log2FC up
  final_up <- merge(global_up_group, 
                    degs_up, by = "Gene_id")
  
  # reordenar columnas
  final_up_clean <- final_up[, c("Chr", "Inicio", "Fin", "log2FC")]
  
  # volver globales
  final_up_glob <- final_up
  final_up_clean_glob <<- final_up_clean
  
  # variable nuevo_dir
  nuevo_dir <- paste(ruta_transcriptoma_modificado, tratamiento, "_circosplot", sep = "")
  
  # crear nuevo directorio
  dir.create(nuevo_dir)
  
  # guardar up
  write.csv(x = final_up_glob,  
            file = paste(nuevo_dir, "/", "circos_up_source.csv", sep = ""),
            row.names = FALSE)
  
  write.csv(x = final_up_clean_glob,  
            file = paste(nuevo_dir, "/", "circos_up.csv", sep = ""),
            row.names = FALSE)
  
  # agregar columna de log2FC down
  final_down <- merge(global_down_group, 
                      degs_down, by = "Gene_id")
  
  # reordenar columnas
  final_down_clean <- final_down[, c("Chr", "Inicio", "Fin", "log2FC")]
  
  # volver globales
  final_down_glob <- final_down
  final_down_clean_glob <<- final_down_clean
  
  # guardar up
  write.csv(x = final_down_glob,  
            file = paste(nuevo_dir, "/", "circos_down_source.csv", sep = ""),
            row.names = FALSE)
  
  write.csv(x = final_down_clean_glob,  
            file = paste(nuevo_dir, "/", "circos_down.csv", sep = ""),
            row.names = FALSE)
  
  ############################# GRAFICAR #################################
  
  # for loop para formatos
  for(i in formatos) {
    
    # crear y guardar los circosplot
    match.fun(i)(paste(paste(nuevo_dir, "/circos_up", sep = ""), ".",i, sep = ""),
                 res = resolucion,
                 width = 3500,
                 height = 3500)
    
    ################ up ###################
    
    ### paramwetros basicos
    circos.par(canvas.xlim = c(-1.2, 1.2), 
               canvas.ylim = c(-1.2, 1.2),
               cell.padding = c(0, 0, 0, 0), 
               start.degree = 90, 
               gap.degree = 4)
    
    # paso 1: inicializar
    circos.genomicInitialize(data = final_up_clean_glob, track.height = 0.15, plotType = NULL)
    
    # paso 2:construir tracks (con colores)
    circos.track(ylim = c(0, 1), factors = final_up_clean_glob$Chr,
                 panel.fun = function(x, y) {
                   chr = CELL_META$sector.index
                   xlim = CELL_META$xlim
                   ylim = CELL_META$ylim
                   
                   # text direction (dd) and adjusmtents (aa)
                   theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
                   dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
                   aa = c(1, 0.5)
                   if(theta < 90 || theta > 270)  aa = c(0, 0.5)
                   
                   # plot chr labels
                   circos.text(x = mean(xlim), y = 1.1, labels = chr, facing = dd, cex = 1, 
                               adj = aa, niceFacing = TRUE)
                   
                   # colores de chr 
                   circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))}, bg.border = NA)
    
    # agregar lineas y puntos 
    circos.genomicTrackPlotRegion(data = final_up_clean_glob, panel.fun = function(x, y, ...) {
      
      # first track, points layout
      circos.genomicLines(x, y, pch = 16, cex = 0.5, col = deg_col_up)
      
      # first track, points layout
      circos.genomicPoints(x, y, pch = 16, cex = 0.5, col = deg_col_up)
      
    })             
  
    
    # limpiar plot
    circos.clear()
    
    dev.off()
  
  } # for loop
  
  ############# down #################
  
  for(i in formatos) {
    
    # crear y guardar los heatmaps
    match.fun(i)(paste(paste(nuevo_dir, "/circos_down", sep = ""), ".",i, sep = ""),
                 res = resolucion,
                 width = 3500,
                 height = 3500)
    
    ### paramwetros basicos
    circos.par(canvas.xlim = c(-1.2, 1.2), 
               canvas.ylim = c(-1.2, 1.2),
               cell.padding = c(0, 0, 0, 0), 
               start.degree = 90, 
               gap.degree = 4)
    
    # paso 1: inicializar
    circos.genomicInitialize(data = final_down_clean_glob, track.height = 0.15, plotType = NULL)
    
    # construir tracks (con colores), configurar primer track
    circos.track(ylim = c(0, 1), factors = final_down_clean_glob$Chr,
                 panel.fun = function(x, y) {
                   chr = CELL_META$sector.index
                   xlim = CELL_META$xlim
                   ylim = CELL_META$ylim
                   
                   # text direction (dd) and adjusmtents (aa)
                   theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
                   dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
                   aa = c(1, 0.5)
                   if(theta < 90 || theta > 270)  aa = c(0, 0.5)
                   
                   # plot chr labels
                   circos.text(x = mean(xlim), y = 1.1, labels = chr, facing = dd, cex = 1, 
                               adj = aa, niceFacing = TRUE)
                   
                   # colores de chr 
                   circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))}, bg.border = NA)
    
    # agregar lineas y puntos
    circos.genomicTrackPlotRegion(data = final_down_clean_glob, panel.fun = function(x, y, ...) {
      
      # second track, lines layout
      circos.genomicLines(x, y, pch = 16, cex = 0.5, col = deg_col_down  )
      
      # second track, points layout
      circos.genomicPoints(x, y, pch = 16, cex = 0.5, col = deg_col_down)
      
    })             
    
    # limpiar plot
    circos.clear()
    
    dev.off()
    
  } # for loop
  
} # function


