test_BUSgibbs <- function(){
    example_Data <- list()
    #batch 1
    example_Data[[1]] <- rbind(matrix(c(1,1,5,5,10,10,
						3,3,7,7,12,12), ncol=6, byrow=TRUE), matrix(c(1,2),nrow=18, ncol=6))

    #batch 2
    batch2_effect <- c(2,2,2,1,1)
    example_Data[[2]] <- rbind(matrix(c(1,1,5,5,10,10,
						3,3,7,7,12,12), ncol=6, byrow=TRUE), matrix(c(1,2),nrow=18, ncol=6)) + batch2_effect

    #batch 3
    batch3_effect <- c(3,2,1,1,2)
    example_Data[[3]] <- rbind(matrix(c(1,1,5,5,10,10,
						3,3,7,7,12,12), ncol=6, byrow=TRUE), matrix(c(1,2),nrow=18, ncol=6)) + batch3_effect

    set.seed(123)
    BUSfits <- BUSgibbs(example_Data, n.subtypes = 3, n.iterations = 100, showIteration = FALSE)

    bat_eff_true0 <- cbind(batch2_effect, batch3_effect)
    bat_eff_true <- rbind(bat_eff_true0, bat_eff_true0, bat_eff_true0, bat_eff_true0)	
    #compare location batch effects
    checkEqualsNumeric(BUSfits$gamma[ ,2:3], bat_eff_true, tol=0.1)	   	
}