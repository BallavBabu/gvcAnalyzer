test_that("BM 2025 Tripartite Trade Identity Holds", {
  # Load internal data
  io <- bm_build_io(bm_toy_Z, bm_toy_Y, bm_toy_VA, bm_toy_X,
                    bm_toy_countries, bm_toy_sectors)

  # Run Trade Decomp
  res <- bm_2025_tripartite_trade(io, 1, 4)

  # Check Identity: E = DAVAX + PF + TS + PB
  calc_E <- res$DAVAX_sr + res$GVC_PF + res$GVC_TS + res$GVC_PB
  expect_equal(res$E_sr, calc_E, tolerance = 1e-8)
})

test_that("BM 2025 Output Identity Holds", {
  io <- bm_build_io(bm_toy_Z, bm_toy_Y, bm_toy_VA, bm_toy_X,
                    bm_toy_countries, bm_toy_sectors)

  res <- bm_2025_output_components(io)

  # Check Identity: X = DomX + TradX + GVC_X
  row1 <- res[1, ]
  calc_X <- row1$DomX + row1$TradX + row1$GVC_X
  expect_equal(row1$X_total, calc_X, tolerance = 1e-8)
})
