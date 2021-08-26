setMethod("ftable",
                      signature(object="metaponeResult"),
                      function(object)
                      {
                              object@mapped.features
                      }
                      )

setMethod("ptable",
                      signature(object="metaponeResult"),
                      function(object)
                      {
                              object@test.result
                      }
                      )
					 