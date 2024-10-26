library(shinytest2)

test_that("{shinytest2} recording: oneptm", {
  app <- AppDriver$new(variant = platform_variant(), name = "oneptm", height = 993, 
      width = 1569)
  app$set_inputs(results_full_screen = FALSE, allow_no_input_binding_ = TRUE)
  app$set_inputs(` Parameter groups` = c("Experimental Design", "Ground Truth Data"))
  app$set_inputs(FracModProt = 1)
  app$click("add_ptm")
  app$set_window_size(width = 1569, height = 993)
  app$set_inputs(overwrite = TRUE)
  app$click("start_sim")
  app$set_window_size(width = 1569, height = 993)
  app$expect_values()
  app$expect_values()
  app$expect_screenshot()
})
