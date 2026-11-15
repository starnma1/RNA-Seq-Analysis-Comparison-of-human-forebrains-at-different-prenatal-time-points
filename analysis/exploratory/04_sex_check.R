all_genes <- gene_counts_annotated_unique
xist_row <- all_genes[all_genes$gene_name == "XIST", ]
rps4y1_row <- all_genes[all_genes$gene_name == "RPS4Y1", ]

stopifnot(nrow(xist_row) == 1, nrow(rps4y1_row) == 1)
stopifnot(identical(names(xist_row), names(rps4y1_row)))

x <- unlist(xist_row[1, 5:ncol(xist_row)], use.names = FALSE)
y <- unlist(rps4y1_row[1, 5:ncol(rps4y1_row)], use.names = FALSE)

samples <- names(xist_row)[5:ncol(xist_row)]
colors <- ifelse(grepl("Female", samples, ignore.case = TRUE), "red",
                 ifelse(grepl("Male", samples, ignore.case = TRUE), "blue", "grey"))

plot(x + 1, y + 1, log = "xy", col = colors, pch = 19,
     xlab = "XIST counts + 1", ylab = "RPS4Y1 counts + 1")

legend("topright",
       legend = c(
         paste0("Female (n=", sum(colors == "red"), ")"),
         paste0("Male (n=", sum(colors == "blue"), ")")
       ),
       col = c("red", "blue"), pch = 19)

