source(testthat::test_path("data-for-tests.R"))

pVar <- mombf::igprior(alpha=0.01, lambda=0.01)
groups <- c(1,2,3,4,5,5,5,6,6,6)

# editable section
pCoef <- mombf::zellnerprior(tau=0.348)
pGroup <- mombf::groupzellnerprior(tau=0.348)

ans_max <- mombf::nlpMarginal(
  theta9_truth_idx, y9, X9, groups=groups, family="normal",
  priorCoef=pCoef, priorVar=pVar, priorGroup=pGroup
)
ans_all <- mombf::nlpMarginal(
  seq_along(theta9_truth), y9, X9, groups=groups, family="normal",
  priorCoef=pCoef, priorVar=pVar, priorGroup=pGroup
)
cat(paste("expected_max=", ans_max, ", expected_all=", ans_all, "\n", sep=""))
