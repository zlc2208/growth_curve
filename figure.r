library(readxl)
library(crayon)
library(ggplot2)
library(dplyr)
library(tidyr)
library(writexl)
library(randomcoloR)
library(gridExtra)
library(ggrepel)
library(colorspace)
library(drc)
library(tibble)
library(ggsignif)
library(rlang)
library(pheatmap)
library(RColorBrewer)
library(effsize)

theme_set(theme_bw())
date <- as.character("20250815")

#----loading growth curve data as od_long----
if (TRUE) {
    if (file.exists(paste(
        "output_", date, "/od_long_", date, ".xlsx",
        sep = ""
    ))) {
        cat(crayon::green("od_long.xlsx exists, loading it...\n"))
        od_long <- read_excel(paste("output_", date, "/od_long_", date, ".xlsx", sep = ""), col_types = "text")
        numeric_col <- sapply(od_long, function(x) all(!is.na(suppressWarnings(as.numeric(x)))))
        od_long[numeric_col] <- lapply(od_long[numeric_col], function(x) as.numeric(x))
        od_long$plasmidid <- paste(od_long$plasmid, od_long$meaning)
    } else {
        cat(crayon::red("od_long.xlsx does not exist\n"))
    }
}
#----loading logistic fit data as gr_lg----
if (TRUE) {
    if (file.exists(paste(
        "output_", date, "/growth_curve_logistic_", date, ".xlsx",
        sep = ""
    ))) {
        cat(crayon::green("growth_curve_logistic.xlsx exists, loading it...\n"))
        gr_lg <- read_excel(paste("output_", date, "/growth_curve_logistic_", date, ".xlsx", sep = ""), col_types = "text")
        numeric_col <- sapply(gr_lg, function(x) all(!is.na(suppressWarnings(as.numeric(x)))))
        gr_lg[numeric_col] <- lapply(gr_lg[numeric_col], function(x) as.numeric(x))
    } else {
        cat(crayon::red("growth_curve_logistic.xlsx does not exist\n"))
    }
}


#----motif mutation growth curve----
if (TRUE) {
    data_control <- filter(od_long, meaning %in% c("positive_control", "negative_control")) %>%
        subset(select = -c(RNA, Protein, meaning))
    data_figure <- od_long[!(is.na(od_long$layout)), ] %>%
        mutate(meaning = factor(meaning, c("WT", "HEXIM1_ILAA", "7SK_m1")))
    if (FALSE) {
        color_figure <- distinctColorPalette(
            length(unique(data_figure$plasmidid)),
            runTsne = TRUE
        )
        names(color_figure) <- unique(data_figure$plasmidid)
        cat(crayon::yellow("color_figure generated successfully.\n"))
        dput(color_figure)
    } else {
        color_figure <- darken(c(
            `pTY128 WT` = "#EEC6E7", `pTY133 WT` = "#CDDE9C", `pTY132 WT` = "#86E86A",
            `pTY134 WT` = "#B127B1", `pTY206 HEXIM1_ILAA` = "#FF93BE", `pTY207 HEXIM1_ILAA` = "#F2DF66",
            `pTY208 HEXIM1_ILAA` = "#E9EEE4", `pTY210 7SK_m1` = "#A0F2C3",
            `pTY209 HEXIM1_ILAA` = "#89ADD9", `pTY211 7SK_m1` = "#FAAD84",
            `pTY212 7SK_m1` = "#92E6F0", `pTY213 7SK_m1` = "#A892F2"
        ), amount = 0.382, fix = TRUE)
    }
    if (FALSE) {
        color_control <- distinctColorPalette(
            length(unique(data_control$plasmidid)),
            runTsne = TRUE
        )
        names(color_control) <- unique(data_control$plasmidid)
        cat(crayon::yellow("color_control generated successfully.\n"))
        dput(color_control)
    } else {
        color_control <- c(
            `pGJJ336 positive_control` = "pink",
            `pTY050 negative_control` = "lightblue"
        )
    }
    color_all <- color_figure
    color_all[names(color_control)] <- color_control
    #----growth curve----
    if (TRUE) {
        p_1 <- ggplot(
            data = data_figure,
            aes(x = time, y = value, group = well, color = plasmidid, alpha = replicate)
        ) +
            facet_grid(RNA + Protein ~ meaning) +
            geom_line() +
            scale_color_manual(values = color_figure) +
            labs(x = "time/h", y = expression(OD[600]))
        ggsave(paste("figure_", date, "/growth_curve_mutant_", date, ".png", sep = ""),
            plot = p_1, width = 11, height = 7.64
        )
        cat(crayon::green("Plot p_1 saved successfully.\n"))
        p_2 <- p_1 + geom_line(
            data = data_control,
            aes(x = time, y = value, group = well, color = plasmidid, alpha = replicate),
            inherit.aes = FALSE,
            linetype = "dashed"
        ) +
            scale_color_manual(values = color_all)
        ggsave(paste("figure_", date, "/growth_curve_mutant_control_", date, ".png", sep = ""),
            plot = p_2, width = 11, height = 7.64
        )
        cat(crayon::green("Plot p_2 saved successfully.\n"))
    }
}
#----motif mutation logistic regression----
if (FALSE) {
    data_figure <- gr_lg[!(is.na(gr_lg$layout)), ] %>%
        mutate(meaning = factor(meaning, c("WT", "HEXIM1_ILAA", "7SK_m1")))
    data_figure$plasmid_type <- paste(data_figure$meaning, data_figure$layout)
    # set color
    if (FALSE) {
        color_figure <- distinctColorPalette(
            length(unique(data_figure$plasmid_type)),
            runTsne = TRUE
        )
        names(color_figure) <- unique(data_figure$plasmid_type)
        cat(crayon::yellow("color_figure generated successfully.\n"))
        dput(color_figure)
    } else {
        color_figure <- darken(c(
            `WT lR3sPN` = "#E2DADE", `WT sR3lPN` = "#C6D691",
            `WT lR3lPN` = "#86E86A", `WT sR3sPN` = "#B755E8",
            `HEXIM1_ILAA sR3sPN` = "#F089AF", `HEXIM1_ILAA sR3lPN` = "#F2DF66",
            `HEXIM1_ILAA lR3sPN` = "#A0F3C8", `7SK_m1 sR3sPN` = "#AFC5F8",
            `HEXIM1_ILAA lR3lPN` = "#D19D6F", `7SK_m1 sR3lPN` = "#92E6F0",
            `7SK_m1 lR3sPN` = "#7493FA", `7SK_m1 lR3lPN` = "#DF86E4"
        ), amount = 0.1, fix = TRUE)
    }
    #----bar_t----
    # odmax
    if (FALSE) {
        odmax_mut_bar <- ggplot(data_figure, aes(x = meaning, y = od_max, fill = plasmid_type)) +
            facet_wrap(~layout, scales = "free_x", nrow = 1) +
            stat_summary(fun = mean, geom = "bar", width = 0.6) + # 绘制条状图
            stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) + # 添加误差线 (标准误)
            scale_fill_manual(values = color_figure) +
            labs(y = expression(OD[600]), title = expression(
                "T-test: " * H[0] * ":" * OD600[max["1"]] * "=" * OD600[max["2"]]
            )) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        for (i in unique(data_figure$layout)) {
            odmax_mut_bar <- odmax_mut_bar +
                geom_signif(
                    data = data_figure[data_figure$layout == i, ],
                    test = "t.test", color = "darkred",
                    comparisons = list(
                        c("HEXIM1_ILAA", "7SK_m1"),
                        c("WT", "HEXIM1_ILAA"),
                        c("WT", "7SK_m1")
                    ),
                    map_signif_level = TRUE, step_increase = 0.1
                )
        }
        ggsave(paste("figure_", date, "/odmax_mut_bar_", date, ".png", sep = ""),
            plot = odmax_mut_bar, width = 10, height = 5
        )
        cat(crayon::green("Plot odmax_single_bar saved successfully.\n"))
        odmax_mut_bar2 <- ggplot(data_figure, aes(x = layout, y = od_max, fill = plasmid_type)) +
            facet_wrap(~meaning, scales = "free_x", nrow = 1) +
            stat_summary(fun = mean, geom = "bar", width = 0.6) + # 绘制条状图
            stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) + # 添加误差线 (标准误)
            scale_fill_manual(values = color_figure) +
            labs(y = expression(OD[600]), title = expression(
                "T-test: " * H[0] * ":" * OD600[max["1"]] * "=" * OD600[max["2"]]
            )) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        for (i in unique(data_figure$meaning)) {
            odmax_mut_bar2 <- odmax_mut_bar2 +
                geom_signif(
                    data = data_figure[data_figure$meaning == i, ],
                    test = "t.test", color = "darkred",
                    comparisons = list(
                        c("lR3lPN", "lR3sPN"),
                        c("sR3lPN", "sR3sPN"),
                        c("lR3sPN", "sR3sPN"),
                        c("lR3lPN", "sR3lPN")
                    ),
                    map_signif_level = TRUE, step_increase = 0.1
                )
        }
        ggsave(paste("figure_", date, "/odmax_mut_bar2_", date, ".png", sep = ""),
            plot = odmax_mut_bar2, width = 10, height = 5
        )
        cat(crayon::green("Plot odmax_single_bar2 saved successfully.\n"))
    }
    # growth_rate
    if (FALSE) {
        growth_rate_mut_bar <- ggplot(data_figure, aes(x = meaning, y = growth_rate, fill = plasmid_type)) +
            facet_wrap(~layout, scales = "free_x", nrow = 1) +
            stat_summary(fun = mean, geom = "bar", width = 0.6) + # 绘制条状图
            stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) + # 添加误差线 (标准误)
            scale_fill_manual(values = color_figure) +
            labs(
                y = expression("Growth rate/(" * OD[600] * "/h)"),
                title = expression(
                    "T-test: " * H[0] * ":" * Growth_rate["1"] * "=" * Growth_rate["2"]
                )
            ) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        for (i in unique(data_figure$layout)) {
            growth_rate_mut_bar <- growth_rate_mut_bar +
                geom_signif(
                    data = data_figure[data_figure$layout == i, ],
                    test = "t.test", color = "darkred",
                    comparisons = list(
                        c("HEXIM1_ILAA", "7SK_m1"),
                        c("WT", "HEXIM1_ILAA"),
                        c("WT", "7SK_m1")
                    ),
                    map_signif_level = TRUE, step_increase = 0.1
                )
        }
        ggsave(paste("figure_", date, "/growth_rate_mut_bar_", date, ".png", sep = ""),
            plot = growth_rate_mut_bar, width = 10, height = 5
        )
        cat(crayon::green("Plot growth_rate_single_bar saved successfully.\n"))
        growth_rate_mut_bar2 <- ggplot(data_figure, aes(x = layout, y = growth_rate, fill = plasmid_type)) +
            facet_wrap(~meaning, scales = "free_x", nrow = 1) +
            stat_summary(fun = mean, geom = "bar", width = 0.6) + # 绘制条状图
            stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) + # 添加误差线 (标准误)
            scale_fill_manual(values = color_figure) +
            labs(
                y = expression("Growth rate/(" * OD[600] * "/h)"),
                title = expression(
                    "T-test: " * H[0] * ":" * Growth_rate["1"] * "=" * Growth_rate["2"]
                )
            ) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        for (i in unique(data_figure$meaning)) {
            growth_rate_mut_bar2 <- growth_rate_mut_bar2 +
                geom_signif(
                    data = data_figure[data_figure$meaning == i, ],
                    test = "t.test", color = "darkred",
                    comparisons = list(
                        c("lR3lPN", "lR3sPN"),
                        c("sR3lPN", "sR3sPN"),
                        c("lR3sPN", "sR3sPN"),
                        c("lR3lPN", "sR3lPN")
                    ),
                    map_signif_level = TRUE, step_increase = 0.1
                )
        }
        ggsave(paste("figure_", date, "/growth_rate_mut_bar2_", date, ".png", sep = ""),
            plot = growth_rate_mut_bar2, width = 10, height = 5
        )
        cat(crayon::green("Plot growth_rate_single_bar2 saved successfully.\n"))
    }
    #---heatmap_MaV----
    if (TRUE) {
        if (TRUE) {
            heatmap_MaV <- function(data_hp = data_hp, variable = variable, rname = rname, cname = cname, type = type, width = 8, height = 8) {
                name_variable <- as_name(ensym(variable))
                rname_row <- as_name(ensym(rname))
                cname_col <- as_name(ensym(cname))
                name_type <- as_name(ensym(type))
                data_MaV <- do.call("rbind", lapply(unique(data_hp[[name_type]]), function(x) {
                    data_hp_sub <- data_hp[data_hp[[name_type]] == x, ]
                    data.frame(
                        rname = unique(data_hp_sub[[rname_row]]),
                        cname = unique(data_hp_sub[[cname_col]]),
                        value_mean = mean(data_hp_sub[[name_variable]]),
                        value_sd = sd(data_hp_sub[[name_variable]])
                    )
                })) %>%
                    mutate(value_CV = value_sd / value_mean * 100) %>%
                    arrange(rname, cname)
                # 平均值矩阵
                mean_mat <- pivot_wider(subset(data_MaV, select = c(rname, cname, value_mean)),
                    names_from = cname, values_from = value_mean
                ) %>%
                    column_to_rownames("rname") %>%
                    as.matrix()
                print(paste("Mean Matrix_", as.character(name_variable), ":", seq = ""))
                print(mean_mat)
                pheatmap(mean_mat,
                    main = paste("Mean of ", as.character(name_variable), seq = ""),
                    cluster_rows = FALSE, cluster_cols = FALSE,
                    color = colorRampPalette(brewer.pal(9, "Blues")[1:7])(255),
                    angle_col = 45, #  热图列名角度
                    display_numbers = TRUE, #  矩阵的数值是否显示在热图上
                    number_format = "%.4f",
                    filename = paste("figure_", date, "/heatmap_mean_",
                        as.character(name_variable), "_", date, ".png",
                        sep = ""
                    ), width = width, height = height
                )
                cat(green(paste("Mean Heatmap_", as.character(name_variable), " finish.\n")))
                # 变异系数矩阵（小10%中100%大）
                CV_mat <- pivot_wider(subset(data_MaV, select = c(rname, cname, value_CV)),
                    names_from = cname,
                    values_from = value_CV,
                    values_fn = mean
                ) %>%
                    column_to_rownames("rname") %>%
                    as.matrix()
                print(paste("CV Matrix_", as.character(name_variable), ":", seq = ""))
                print(CV_mat)
                breaks <- c(seq(0, 34, length.out = 255), 35, seq(36, ifelse(max(data_MaV$value_CV) > 100, max(data_MaV$value_CV), 100), length.out = 255))
                colors <- c(colorRampPalette(rev(brewer.pal(9, "Blues")[1:6]))(255), "white", colorRampPalette(brewer.pal(9, "Reds")[1:7])(255))
                pheatmap(CV_mat,
                    main = paste("CV% of ", as.character(name_variable), seq = ""),
                    cluster_rows = FALSE, cluster_cols = FALSE,
                    breaks = breaks,
                    color = colors,
                    angle_col = 45, #  热图列名角度
                    display_numbers = TRUE, #  矩阵的数值是否显示在热图上
                    number_format = "%.2f",
                    filename = paste("figure_", date,
                        "/heatmap_CV_", as.character(name_variable), "_", date, ".png",
                        sep = ""
                    ), width = width, height = height
                )
                cat(green(paste("CV Heatmap_", as.character(name_variable), " finish.\n")))
            }
            save(heatmap_MaV, file = "d:/app/r/function/heatmap_MaV.RDate")
        } else {
            load("d:/app/r/function/heatmap_MaV.RDate")
        }
        heatmap_MaV(data_hp = data_figure, variable = od_max, rname = layout, cname = meaning, type = plasmid_type, width = 5, height = 5)
        heatmap_MaV(data_hp = data_figure, variable = growth_rate, rname = layout, cname = meaning, type = plasmid_type, width = 5, height = 5)
    }
    #----heatmap_ttest----
    if (FALSE) {
        data_test <- subset(data_figure, select = c(plasmid_type, od_max, growth_rate, meaning, layout))
        types <- unique(data_test$plasmid_type)
        diff_pairs <- combn(types, 2, simplify = FALSE)
        diff_pairs_rev <- lapply(diff_pairs, function(pair) rev(pair))
        self_pairs <- lapply(types, function(x) c(x, x))
        all_pairs <- c(diff_pairs, diff_pairs_rev, self_pairs)
        #---- odmax----
        t_odmax <- do.call(rbind, lapply(all_pairs, function(pair) {
            test <- t.test(
                data_test$od_max[data_test$plasmid_type == pair[1]],
                data_test$od_max[data_test$plasmid_type == pair[2]],
                alternative = "greater"
            )
            cd <- cohen.d(
                data_test$od_max[data_test$plasmid_type == pair[1]],
                data_test$od_max[data_test$plasmid_type == pair[2]]
            )
            data.frame(
                layout1 = unique(data_test$layout[data_test$plasmid_type == pair[1]]),
                layout2 = unique(data_test$layout[data_test$plasmid_type == pair[2]]),
                meaning1 = unique(data_test$meaning[data_test$plasmid_type == pair[1]]),
                meaning2 = unique(data_test$meaning[data_test$plasmid_type == pair[2]]),
                p_value = test$p.value,
                cohen_d = cd$estimate
            )
        })) %>%
            as.data.frame() %>%
            mutate(
                p_value = as.numeric(p_value),
                cohen_d = as.numeric(cohen_d)
                # DHFR_position1 = factor(DHFR_position1, level = (unique(DHFR_position1))),
                # DHFR_position2 = factor(DHFR_position2, level = rev(unique(DHFR_position1)))
            )
        # set color bar
        color_pvalue <- colorRampPalette(rev(brewer.pal(9, "YlGnBu")[1:6]))(255)
        color_cohend <- colorRampPalette((brewer.pal(11, "RdBu")[2:10]))(255)
        p_t_odmax_mut <- ggplot(data = t_odmax, aes(x = layout2, y = layout1)) +
            facet_grid(meaning1 ~ meaning2, scales = "fixed") +
            coord_fixed(ratio = 1) +
            geom_tile(aes(fill = p_value)) +
            # geom_text(aes(label = format(p_value, scientific = TRUE, digits = 3))) +
            scale_fill_gradientn(colours = color_pvalue) +
            labs(title = expression("T-test: " * H[0] * ":" * OD600[max[layout1]] - OD600[max[layout2]] <= 0)) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        ggsave(
            plot = p_t_odmax_mut,
            filename = paste("figure_", date, "/t_test_odmax_mut.png", sep = ""),
            width = 11, height = 9
        )
        print("t_test_odmax_mut finish.")
        p_cd_odmax_mut <- ggplot(data = t_odmax, aes(x = layout2, y = layout1)) +
            facet_grid(meaning1 ~ meaning2, scales = "fixed") +
            coord_fixed(ratio = 1) +
            geom_tile(aes(fill = cohen_d)) +
            # geom_text(aes(label = round(cohen_d, 2))) +
            scale_fill_gradientn(colours = color_cohend) +
            labs(title = expression("Cohen's d:" * OD600[max[layout1]] - OD600[max[layout2]])) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        ggsave(
            plot = p_cd_odmax_mut,
            filename = paste("figure_", date, "/cohen_d_odmax_mut.png", sep = ""),
            width = 11, height = 9
        )
        print("cohen_d_odmax_mut finish.")
        #---- growth_rate----
        t_gr <- do.call(rbind, lapply(all_pairs, function(pair) {
            test <- t.test(
                data_test$growth_rate[data_test$plasmid_type == pair[1]],
                data_test$growth_rate[data_test$plasmid_type == pair[2]],
                alternative = "greater"
            )
            cd <- cohen.d(
                data_test$growth_rate[data_test$plasmid_type == pair[1]],
                data_test$growth_rate[data_test$plasmid_type == pair[2]]
            )
            data.frame(
                layout1 = unique(data_test$layout[data_test$plasmid_type == pair[1]]),
                layout2 = unique(data_test$layout[data_test$plasmid_type == pair[2]]),
                meaning1 = unique(data_test$meaning[data_test$plasmid_type == pair[1]]),
                meaning2 = unique(data_test$meaning[data_test$plasmid_type == pair[2]]),
                p_value = test$p.value,
                cohen_d = cd$estimate
            )
        })) %>%
            as.data.frame() %>%
            mutate(
                p_value = as.numeric(p_value),
                cohen_d = as.numeric(cohen_d)
                # DHFR_position1 = factor(DHFR_position1, levels = (unique(DHFR_position1))),
                # DHFR_position2 = factor(DHFR_position2, levels = rev(unique(DHFR_position1)))
            )
        p_t_gr_mut <- ggplot(data = t_gr, aes(x = layout2, y = layout1)) +
            facet_grid(meaning1 ~ meaning2, scales = "fixed") +
            coord_fixed(ratio = 1) +
            geom_tile(aes(fill = p_value)) +
            # geom_text(aes(label = format(p_value, scientific = TRUE, digits = 3))) +
            scale_fill_gradientn(colours = color_pvalue) +
            labs(title = expression("T-test: " * H[0] * ":" * Growth_rate[layout1] - Growth_rate[layout2] <= 0)) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        ggsave(
            plot = p_t_gr_mut,
            filename = paste("figure_", date, "/t_test_gr_mut.png", sep = ""),
            width = 11, height = 9
        )
        print("t_test_gr_mut finish.")
        p_cd_gr_mut <- ggplot(data = t_gr, aes(x = layout2, y = layout1)) +
            facet_grid(meaning1 ~ meaning2, scales = "fixed") +
            coord_fixed(ratio = 1) +
            geom_tile(aes(fill = cohen_d)) +
            # geom_text(aes(label = round(cohen_d, 2))) +
            scale_fill_gradientn(colours = color_cohend) +
            labs(title = expression("Cohen's d:" * Growth_rate[layout1] - Growth_rate[layout2])) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        ggsave(
            plot = p_cd_gr_mut,
            filename = paste("figure_", date, "/cohen_d_gr_mut.png", sep = ""),
            width = 11, height = 9
        )
        print("cohen_d_gr_mut finish.")
    }
}


#----single positive growth curve----
if (FALSE) {
    data_figure <- od_long[od_long$meaning %in% c("double positive", "single protein", "single RNA", "double negative"), ] %>%
        mutate(meaning = factor(meaning, levels = c("double positive", "single protein", "single RNA", "double negative")))
    data_control <- filter(od_long, meaning %in% c("positive_control", "negative_control")) %>%
        subset(select = -c(meaning))
    # set color
    if (FALSE) {
        color_figure <- distinctColorPalette(
            length(unique(data_figure$plasmidid)),
            runTsne = TRUE
        )
        names(color_figure) <- unique(data_figure$plasmidid)
        cat(crayon::yellow("color_figure generated successfully.\n"))
        dput(color_figure)
    } else {
        color_figure <- darken(c(
            `pTY128 double positive` = "#C876D5", `pTY135 single protein` = "#95C15B",
            `pTY136 single RNA` = "#9EBAC0", `pTY121 double negative` = "#FFC4AC"
        ), amount = 0.382, fix = TRUE)
    }
    if (FALSE) {
        color_control <- distinctColorPalette(
            length(unique(data_control$plasmidid)),
            runTsne = TRUE
        )
        names(color_control) <- unique(data_control$plasmidid)
        cat(crayon::yellow("color_control generated successfully.\n"))
        dput(color_control)
    } else {
        color_control <- c(
            `pGJJ336 positive_control` = "pink",
            `pTY050 negative_control` = "lightblue"
        )
    }
    color_all <- color_figure
    color_all[names(color_control)] <- color_control
    # growth curve
    if (FALSE) {
        p_single <- ggplot(
            data_figure,
            aes(x = time, y = value, group = well, color = plasmidid, alpha = replicate)
        ) +
            facet_wrap(~meaning, nrow = 1) +
            geom_line() +
            scale_color_manual(values = color_figure) +
            labs(x = "time/h", y = expression(OD[600]))
        ggsave(paste(
            "figure_", date, "/growth_curve_single_", date, ".png",
            sep = ""
        ), plot = p_single, width = 12, height = 4)
        cat(crayon::green("Plot p_single saved successfully.\n"))
        p_single_control <- p_single +
            geom_line(
                data = data_control,
                aes(x = time, y = value, group = well, color = plasmidid, alpha = replicate),
                inherit.aes = FALSE,
                linetype = "dashed"
            ) +
            scale_color_manual(values = color_all)
        ggsave(paste(
            "figure_", date, "/growth_curve_single_control_", date, ".png",
            sep = ""
        ), plot = p_single_control, width = 12, height = 4)
        cat(crayon::green("Plot p_single_control saved successfully.\n"))
    }
}
#----single positive logistic regression----
if (FALSE) {
    data_figure <- gr_lg[gr_lg$meaning %in% c("double positive", "single protein", "single RNA", "double negative"), ] %>%
        mutate(meaning = factor(meaning, levels = c("double positive", "single protein", "single RNA", "double negative")))
    # set color
    if (FALSE) {
        color_figure <- distinctColorPalette(
            length(unique(data_figure$meaning)),
            runTsne = TRUE
        )
        names(color_figure) <- unique(data_figure$meaning)
        cat(crayon::yellow("color_figure generated successfully.\n"))
        dput(color_figure)
    } else {
        color_figure <- darken(c(
            `double positive` = "#C876D5", `single protein` = "#95C15B",
            `single RNA` = "#9EBAC0", `double negative` = "#FFC4AC"
        ), amount = 0.1, fix = TRUE)
    }
    #----bar_t----
    if (FALSE) {
        # odmax
        odmax_single_bar <- ggplot(data_figure, aes(x = meaning, y = od_max, fill = meaning)) +
            stat_summary(fun = mean, geom = "bar", width = 0.6) + # 绘制条状图
            stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) + # 添加误差线 (标准误)
            scale_fill_manual(values = color_figure) +
            geom_signif(
                test = "t.test", color = "darkred",
                comparisons = list(
                    c("double positive", "single protein"),
                    c("double positive", "single RNA"),
                    c("double positive", "double negative")
                ),
                map_signif_level = TRUE, step_increase = 0.1
            ) +
            labs(y = expression(OD[600]), title = expression(
                "T-test: " * H[0] * ":" * OD600[max["1"]] * "=" * OD600[max["2"]]
            ))
        ggsave(paste("figure_", date, "/odmax_single_bar_", date, ".png", sep = ""),
            plot = odmax_single_bar, width = 10, height = 5
        )
        cat(crayon::green("Plot odmax_single_bar saved successfully.\n"))
        # growth_rate
        growth_rate_single_bar <- ggplot(data_figure, aes(x = meaning, y = growth_rate, fill = meaning)) +
            stat_summary(fun = mean, geom = "bar", width = 0.6) + # 绘制条状图
            stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) + # 添加误差线 (标准误)
            scale_fill_manual(values = color_figure) +
            geom_signif(
                test = "t.test", color = "darkred",
                comparisons = list(
                    c("double positive", "single protein"),
                    c("double positive", "single RNA"),
                    c("double positive", "double negative")
                ),
                map_signif_level = TRUE, step_increase = 0.1
            ) +
            labs(
                y = expression("Growth rate/(" * OD[600] * "/h)"),
                title = expression(
                    "T-test: " * H[0] * ":" * Growth_rate["1"] * "=" * Growth_rate["2"]
                )
            )
        ggsave(paste("figure_", date, "/growth_rate_single_bar_", date, ".png", sep = ""),
            plot = growth_rate_single_bar, width = 10, height = 5
        )
        cat(crayon::green("Plot growth_rate_single_bar saved successfully.\n"))
    }
    #----heatmap_ttest----
    if (TRUE) {
        data_test <- subset(data_figure, select = c(od_max, growth_rate, meaning))
        types <- unique(data_test$meaning)
        diff_pairs <- combn(types, 2, simplify = FALSE)
        diff_pairs_rev <- lapply(diff_pairs, function(pair) rev(pair))
        self_pairs <- lapply(types, function(x) c(x, x))
        all_pairs <- c(diff_pairs, diff_pairs_rev, self_pairs)
        #---- odmax----
        t_odmax <- do.call(rbind, lapply(all_pairs, function(pair) {
            test <- t.test(
                data_test$od_max[data_test$meaning == pair[1]],
                data_test$od_max[data_test$meaning == pair[2]],
                alternative = "greater"
            )
            cd <- cohen.d(
                data_test$od_max[data_test$meaning == pair[1]],
                data_test$od_max[data_test$meaning == pair[2]]
            )
            data.frame(
                meaning1 = unique(data_test$meaning[data_test$meaning == pair[1]]),
                meaning2 = unique(data_test$meaning[data_test$meaning == pair[2]]),
                p_value = test$p.value,
                cohen_d = cd$estimate
            )
        })) %>%
            as.data.frame() %>%
            mutate(
                p_value = as.numeric(p_value),
                cohen_d = as.numeric(cohen_d)
            )
        # set color bar
        color_pvalue <- colorRampPalette(rev(brewer.pal(9, "YlGnBu")[1:6]))(255)
        color_cohend <- colorRampPalette((brewer.pal(11, "RdBu")[2:10]))(255)
        p_t_odmax_single <- ggplot(data = t_odmax, aes(x = meaning2, y = meaning1)) +
            coord_fixed(ratio = 1) +
            geom_tile(aes(fill = p_value)) +
            geom_text(aes(label = format(p_value, scientific = TRUE, digits = 3))) +
            scale_fill_gradientn(colours = color_pvalue) +
            labs(title = expression("T-test: " * H[0] * ":" * OD600[max[meaning1]] - OD600[max[meaning2]] <= 0)) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        ggsave(
            plot = p_t_odmax_single,
            filename = paste("figure_", date, "/t_test_odmax_single.png", sep = ""),
            width = 7, height = 6
        )
        print("t_test_odmax_single finish.")
        p_cd_odmax_single <- ggplot(data = t_odmax, aes(x = meaning2, y = meaning1)) +
            coord_fixed(ratio = 1) +
            geom_tile(aes(fill = cohen_d)) +
            geom_text(aes(label = round(cohen_d, 2))) +
            scale_fill_gradientn(colours = color_cohend) +
            labs(title = expression("Cohen's d:" * OD600[max[meaning1]] - OD600[max[meaning2]])) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        ggsave(
            plot = p_cd_odmax_single,
            filename = paste("figure_", date, "/cohen_d_odmax_single.png", sep = ""),
            width = 7, height = 6
        )
        print("cohen_d_odmax_single finish.")
        #---- growth_rate----
        t_gr <- do.call(rbind, lapply(all_pairs, function(pair) {
            test <- t.test(
                data_test$growth_rate[data_test$meaning == pair[1]],
                data_test$growth_rate[data_test$meaning == pair[2]],
                alternative = "greater"
            )
            cd <- cohen.d(
                data_test$growth_rate[data_test$meaning == pair[1]],
                data_test$growth_rate[data_test$meaning == pair[2]]
            )
            data.frame(
                meaning1 = unique(data_test$meaning[data_test$meaning == pair[1]]),
                meaning2 = unique(data_test$meaning[data_test$meaning == pair[2]]),
                p_value = test$p.value,
                cohen_d = cd$estimate
            )
        })) %>%
            as.data.frame() %>%
            mutate(
                p_value = as.numeric(p_value),
                cohen_d = as.numeric(cohen_d)
            )
        p_t_gr_single <- ggplot(data = t_gr, aes(x = meaning2, y = meaning1)) +
            coord_fixed(ratio = 1) +
            geom_tile(aes(fill = p_value)) +
            geom_text(aes(label = format(p_value, scientific = TRUE, digits = 3))) +
            scale_fill_gradientn(colours = color_pvalue) +
            labs(title = expression("T-test: " * H[0] * ":" * Growth_rate[meaning1] - Growth_rate[meaning2] <= 0)) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        ggsave(
            plot = p_t_gr_single,
            filename = paste("figure_", date, "/t_test_gr_single.png", sep = ""),
            width = 7, height = 6
        )
        print("t_test_gr_single finish.")
        p_cd_gr_single <- ggplot(data = t_gr, aes(x = meaning2, y = meaning1)) +
            coord_fixed(ratio = 1) +
            geom_tile(aes(fill = cohen_d)) +
            geom_text(aes(label = round(cohen_d, 2))) +
            scale_fill_gradientn(colours = color_cohend) +
            labs(title = expression("Cohen's d:" * Growth_rate[meaning1] - Growth_rate[meaning2])) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        ggsave(
            plot = p_cd_gr_single,
            filename = paste("figure_", date, "/cohen_d_gr_single.png", sep = ""),
            width = 7, height = 6
        )
        print("cohen_d_gr_single finish.")
    }
}
