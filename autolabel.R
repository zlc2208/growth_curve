library(readxl)
library(dplyr)
library(writexl)
library(tidyr)
library(openxlsx)
library(crayon)

date <- "20250815"
file_name <- paste0("label_", date, ".xlsx")

if (file.exists(file_name)) {
    wb <- loadWorkbook(file_name)
    sheet_list <- excel_sheets(file_name)
    # 变换布局图
    if (TRUE) {
        data_map_origin <- read_excel(file_name, sheet = "map_origin") %>%
            mutate(across(everything(), as.character))
        data_map <- data_map_origin
        data_map[, seq(2, 12, by = 2)] <- data_map_origin[, seq(2, 7, by = 1)]
        data_map[, seq(3, 13, by = 2)] <- data_map_origin[, seq(8, 13, by = 1)]
        if ("map" %in% sheet_list) {
            cat(crayon::red(("Sheet map already exists. Removing it.\n")))
            removeWorksheet(wb, "map") # 删除旧工作表
        } else {
            cat(crayon::green(("Sheet map does not exist. Creating a new one.\n")))
        }
        # 创建新的工作簿
        addWorksheet(wb, "map")
        writeData(wb, sheet = "map", x = data_map) # 正确写入数据
        saveWorkbook(wb, file = file_name, overwrite = TRUE)
        cat(crayon::green(("map has been saved.\n")))
    }
    # auto_label based on map
    if (TRUE) {
        data_map <- read_excel(file_name, sheet = "map") %>%
            mutate(across(everything(), as.character))
        data_property <- read_excel(file_name, sheet = "property") %>%
            mutate(across(everything(), as.character))
        label_map <- pivot_longer(
            data_map,
            cols = -c("row"), names_to = "col", values_to = "plasmid"
        )
        label_map$well <- paste(label_map$row, label_map$col, sep = "")
        label_map <- label_map %>% filter(!is.na(plasmid))
        index_nublank <- label_map$plasmid != "blank"
        label_map$plasmid[index_nublank] <- paste("pTY", label_map$plasmid[index_nublank], sep = "")
        auto_label <- merge(label_map, data_property, by = "plasmid")
        auto_label$plasmidid <- paste(auto_label$plasmid, auto_label$meaning, sep = "_")
        blank_row <- subset(label_map, subset = (plasmid == "blank"))
        auto_label <- bind_rows(auto_label, blank_row)
        auto_label <- subset(auto_label, select = -c(row, col))
        auto_label <- auto_label %>%
            group_by(plasmidid) %>%
            mutate(replicate = as.character(row_number())) %>%
            ungroup()
        if ("auto_label" %in% sheet_list) {
            cat(crayon::red(("Sheet auto_label already exists. Removing it.\n")))
            removeWorksheet(wb, "auto_label") # 删除旧工作表
        } else {
            cat(crayon::green(("Sheet auto_label does not exist. Creating a new one.\n")))
        }
        # 创建新的工作簿
        addWorksheet(wb, "auto_label")
        writeData(wb, sheet = "auto_label", x = auto_label) # 正确写入数据
        saveWorkbook(wb, file = file_name, overwrite = TRUE)
        cat(crayon::green(("auto_Label has been saved.\n")))
    }
} else {
    cat(crayon::red(("label.xlsx does not exist. Please check the file path.\n")))
}
