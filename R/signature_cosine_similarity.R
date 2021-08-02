#Calculating the cosine similarity between a vector of mutations and a mutational signature

signature_cosine_similarity <- function(mutations.vector, signature.matrix){
  #create a table of probabilities for the vector of mutations
  tab <- table(mutations.vector)/sum(table(mutations.vector))
  #transform table to a data frame and assign column names to mutations and frequencies
  mutations.df <- data.frame(tab)
  colnames(mutations.df) <- c("mutations", "frequencies")

  #defines a function that joins two data frames and replaces NA with 0
  left_join_NA <- function(x, y, by) {
    left_join(x = x, y = y, by = by) %>%
      mutate_each(funs(replace(., which(is.na(.)), 0)))
  }

  #creates a vector of the 96 mutation types
  mutations <- get_mutation_types()
  #converts vector of 96 mutation types to a data frame
  mutation.types.df <- data.frame(mutations)
  #joins together sample probabilities with data frame of mutation types so that all 96 mutation types are included
  complete.mutations.df <- left_join_NA(mutation.types.df, mutations.df, by = "mutations")

  #converts data frame of mutation types and frequencies to a data frame and assigns rownames
  mutations.matrix <- as.matrix(complete.mutations.df$frequencies)
  rownames(mutations.matrix) = mutations

  #calculates and returns cosine similarity to the matrix for the signature of interest
  cosine.similarity.matrix <- cos_sim_matrix(signature.matrix, mutations.matrix)
  return(cosine.similarity.matrix[,1])
}
