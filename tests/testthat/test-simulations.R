library(assertthat)
context("simulations")
rng <- sort(rgamma(2, shape = 1, scale = 2)) #Range of eigen values
num.vars <- sample(c(5:8), 1) #Number of variables
method <- sample(c("ordered", "ic-dag", "melancon"), 1) #Graph sim method
net <- simGaussianNet(num.vars, rng, method)

#Data frame tests
testSimData <- function(out.data){
  expect_true(is.data.frame(out.data))
  expect_true(noNA(out.data))
  out.data %>% unlist %>% is.numeric %>% all %>% expect_true
  out.data %>% unlist %>% is.finite %>% all %>% expect_true
}

test_that("Output matrices are symetric positive definate.", {
  net %>% bn.net %>% moral %>% 
    getPrecisionBasis %>% simSparsePrecision(rng) %>% eigen %$%
    values %>% min %>%
    expect_more_than(0)
  net %>% bn.net %>% simNetCovariance(rng) %>% eigen %$% values %>% min %>%
    expect_more_than(0)
})
test_that("Output matrices have desired eigen range.", {
  #Output of simSparsePrecision
  basis <- bn.net(net) %>% getPrecisionBasis 
  eigen.range <- simSparsePrecision(basis, rng) %>% eigen %$% values %>% range
  expect_more_than(eigen.range[1], rng[1])
  expect_less_than(eigen.range[2], rng[2])
  #Output of simNetCovariance
  eigen.range <- simNetCovariance(bn.net(net), rng) %>% eigen %$% values %>% range 
  expect_more_than(eigen.range[1], rng[1])
  expect_less_than(eigen.range[2], rng[2])
})

test_that("Output of convert2BinaryMat has named rows and columns", {
  d.names <- clusterGeneration::genPositiveDefMat(10, "unifcorrmat")$Sigma %>%
    convert2BinaryMat %>%
    dimnames
  expect_true(length(d.names) > 0)
})

test_that("Output from Bayesian network simulation can simulate data.", {
  out.data <- net %>%
    rbn(1000) %>%
    testSimData
  
  out.data <- clusterGeneration::genPositiveDefMat(10, "unifcorrmat")$Sigma %>%
    simBNFromCovMat %>%
    rbn(1000) %>%
    testSimData
})

test_that("Moral graph modeling works as expected, since simulation is based on creation of undirected Markov nets (moral graphs) from 
          matrices", {
            #Get a precision basis from the input network via the moral graph, then try to reconstruct
            #the moral graph.
            moral.net.original <- moral(bn.net(net))
            net.out <- moral.net.original %>% #Start with original 
              getPrecisionBasis %>% #Get the basis of the precision
              `diag<-`(0) %>% #Make this aan adjacency matrix by converting diagonals to 0
              adjMat2Moral #Convert this back to a bn graph object
            #This should be a moral network, so...
            expect_equal(net.out %>% directed.arcs %>% nrow, 0) # ... confirm no directed arcs
            expect_more_than(net.out %>% undirected.arcs %>% nrow, 0) #confirm some undirected arcs
            expect_true(all.equal(moral.net.original, net.out)) # input and output are equivalent
            
            #Sim a new graph from the moral graph of another.  Make sure the moral graph of the new one
            #matches the first.
            moral.net.original %>% #Start with original moral net
              getPrecisionBasis %>% #get a basis of precision matrix
              simSparsePrecision(rng) %>% #simulate a valid precision matrix
              solve %>% #take inverse to get valid covariance structure
              simDAGFromCovMat %>% #Simulate a new DAG consistent with covariance structure
              moral %>% # Get the moral network of the new DAG
              all.equal(moral.net.original) %>% #Compare to original moral net
              expect_true
            })

test_that("simDAGFromCovMat should return NULL for a net that is too complex", {
  randomNet(num.vars = 40, method = "melancon") %>%
    simNetCovariance(c(1.5, 5)) %>%
    simDAGFromCovMat %>%
    is.null %>% 
    expect_true
})
