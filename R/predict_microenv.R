#' @importFrom randomForest randomForest
NULL
#' Predict microenvironment cell type
#'
#' This function uses a Random Forest model trained on synthetic
#' gene expression and pathway features to predict the microenvironment
#' cell type (Cancer / T_Cell / Fibroblast).
#'
#' @param newdata A data frame containing at least the following
#'   numeric columns:
#'   \code{Gene_E_Housekeeping}, \code{Gene_A_Oncogene},
#'   \code{Gene_B_Immune}, \code{Gene_C_Stromal},
#'   \code{Gene_D_Therapy}, \code{Pathway_Score_Inflam}.
#'
#' @return A factor vector of predicted \code{Cell_Type}.
#' @examples
#' \dontrun{
#' # Example (assuming you have a data frame 'df' with the required columns):
#' preds <- predict_cell_type(df)
#' table(preds)
#' }
#' @export
predict_cell_type <- function(newdata) {
  required_cols <- c(
    "Gene_E_Housekeeping",
    "Gene_A_Oncogene",
    "Gene_B_Immune",
    "Gene_C_Stromal",
    "Gene_D_Therapy",
    "Pathway_Score_Inflam"
  )

  missing_cols <- setdiff(required_cols, names(newdata))
  if (length(missing_cols) > 0) {
    stop(
      "The following required columns are missing from 'newdata': ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # 保证是 data.frame
  newdata <- as.data.frame(newdata)

  # 只取需要的列（并按训练时顺序排列）
  x <- newdata[, required_cols]

  # 使用内部的 rf_celltype 对象进行预测
  preds <- stats::predict(rf_celltype, newdata = x, type = "class")

  return(preds)
}



#' Predict disease status (Tumor vs Healthy)
#'
#' This function uses a Random Forest model trained on gene expression
#' and engineered pathway-like scores to predict disease status
#' (Tumor vs Healthy_Control).
#'
#' Internally, the function computes the following scores for each row:
#' \itemize{
#'   \item Tumor_score   = Gene_A_Oncogene + Gene_D_Therapy
#'   \item Immune_score  = Gene_B_Immune    + Pathway_Score_Inflam
#'   \item Stromal_score = Gene_C_Stromal
#' }
#'
#' @param newdata A data frame containing at least the following
#'   numeric columns:
#'   \code{Gene_E_Housekeeping}, \code{Gene_A_Oncogene},
#'   \code{Gene_B_Immune}, \code{Gene_C_Stromal},
#'   \code{Gene_D_Therapy}, \code{Pathway_Score_Inflam}.
#' @param type Character string, either \code{"class"} (default)
#'   for predicted labels, or \code{"prob"} for class probabilities.
#'
#' @return If \code{type = "class"}, a factor vector of predicted
#'   \code{Disease_Status}. If \code{type = "prob"}, a data frame
#'   with predicted class probabilities.
#'
#' @examples
#' \dontrun{
#' preds <- predict_disease_status(df)
#' head(preds)
#'
#' probs <- predict_disease_status(df, type = "prob")
#' head(probs)
#' }
#' @export
predict_disease_status <- function(newdata, type = c("class", "prob")) {
  type <- match.arg(type)

  required_cols <- c(
    "Gene_E_Housekeeping",
    "Gene_A_Oncogene",
    "Gene_B_Immune",
    "Gene_C_Stromal",
    "Gene_D_Therapy",
    "Pathway_Score_Inflam"
  )

  missing_cols <- setdiff(required_cols, names(newdata))
  if (length(missing_cols) > 0) {
    stop(
      "The following required columns are missing from 'newdata': ",
      paste(missing_cols, collapse = ", ")
    )
  }

  newdata <- as.data.frame(newdata)

  # 先保留必需列
  newdata <- newdata[, required_cols]

  # 在函数内部复现你们训练脚本里的特征工程
  newdata$Tumor_score   <- newdata$Gene_A_Oncogene + newdata$Gene_D_Therapy
  newdata$Immune_score  <- newdata$Gene_B_Immune   + newdata$Pathway_Score_Inflam
  newdata$Stromal_score <- newdata$Gene_C_Stromal

  # 用内部的 rf_disease_A_score 模型预测
  if (type == "class") {
    preds <- stats::predict(rf_disease_A_score, newdata = newdata, type = "class")
    return(preds)
  } else {
    prob  <- stats::predict(rf_disease_A_score, newdata = newdata, type = "prob")
    return(as.data.frame(prob))
  }
}
#' Predict microenvironment outcomes (labels only)
#'
#' A unified prediction interface for the microenvRF package.
#' For classification tasks, this function returns class labels (not probabilities),
#' as required by the BIO215 coursework.
#'
#' @param newdata A data.frame containing required feature columns.
#'   Required columns for both tasks:
#'   \code{Gene_E_Housekeeping}, \code{Gene_A_Oncogene}, \code{Gene_B_Immune},
#'   \code{Gene_C_Stromal}, \code{Gene_D_Therapy}, \code{Pathway_Score_Inflam}.
#' @param task Character, either \code{"cell_type"} or \code{"disease_status"}.
#'
#' @return A factor vector of predicted class labels.
#' @examples
#' \dontrun{
#' newdata <- data.frame(
#'   Gene_E_Housekeeping = rnorm(5),
#'   Gene_A_Oncogene     = rnorm(5),
#'   Gene_B_Immune       = rnorm(5),
#'   Gene_C_Stromal      = rnorm(5),
#'   Gene_D_Therapy      = rnorm(5),
#'   Pathway_Score_Inflam= rnorm(5)
#' )
#' predict_microenv(newdata, task = "cell_type")
#' predict_microenv(newdata, task = "disease_status")
#' }
#' @export
predict_microenv <- function(newdata, task = c("cell_type", "disease_status")) {
  task <- match.arg(task)

  # 统一检查必需列
  required_cols <- c(
    "Gene_E_Housekeeping",
    "Gene_A_Oncogene",
    "Gene_B_Immune",
    "Gene_C_Stromal",
    "Gene_D_Therapy",
    "Pathway_Score_Inflam"
  )

  missing_cols <- setdiff(required_cols, names(newdata))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in 'newdata': ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # 两个任务都只返回标签（class）
  if (task == "cell_type") {
    return(predict_cell_type(newdata))
  } else {
    # disease_status 默认返回标签
    return(predict_disease_status(newdata, type = "class"))
  }
}
