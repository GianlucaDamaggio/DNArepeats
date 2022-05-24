#' @title DNArepeats
#'
#' @description Instability & Expansion Index from sized DNA repeats; R package for calculating instability/expansion index from DNA repeated region already sized. Works on time-series.
#'
#
#' @param dataframe, threshold
#'
#' @return dataframe
#'
#' @examples
#'
#' @export


expansion_index <- function(dataframe, threshold)
{
  temp_df_to_build<- data.frame(matrix(ncol = 10, nrow = 0))

  for( s in unique(dataframe$sample)){

    for( d in unique(dataframe$day)){

      df_count = dataframe %>% filter(!is.na(count)) %>% filter(day == d , count  >= count[which(n == max(n))])

      df_thresh = df_count %>% filter( n > max(df_count$n) * threshold /100 )

      df_norm_peaks = df_thresh %>% mutate(norm_peaks= n/sum(df_thresh$n)) %>% data.frame

      df_steps = df_norm_peaks %>% mutate(steps = count - (min(count)))

      df_norm_steps = df_steps %>% mutate(norm_steps = norm_peaks * steps)

      df_expansion_index = df_norm_steps %>% mutate(expansion_index = sum(norm_steps))

      df_expansion_index$threshold = threshold

      temp_df_to_build = rbind(temp_df_to_build,df_expansion_index)
    }
  }

  dataframe = temp_df_to_build

  return(dataframe)

}





instability_index <- function(dataframe, threshold)
{

  temp_df_to_build<- data.frame(matrix(ncol = 10, nrow = 0))

  for( s in unique(dataframe$sample)){

    starting_count = as.integer(dataframe %>% filter(count > 0 , day == 0 , sample == s ) %>% select(count, n) %>% filter(count == count[which(n == max(n))]) %>% select(count))

    for( d in unique(dataframe$day)){

      df_count =  dataframe %>% filter(!is.na(count), day == d)

      df_count$steps= -(which(df_count$count == (starting_count - 1))):(nrow(df_count)-(which(df_count$count == starting_count)))

      df_thresh = df_count %>% filter( n > max(df_count$n) * threshold / 100 )

      df_norm_peaks = df_thresh %>% mutate(norm_peaks = n/sum(df_thresh$n)) %>% data.frame

      df_norm_steps = df_norm_peaks %>% mutate(norm_steps = norm_peaks * steps)

      df_instability_index = df_norm_steps %>% mutate(instability_index = sum(norm_steps))

      df_instability_index$threshold = threshold

      temp_df_to_build = rbind(temp_df_to_build,df_instability_index)
    }
  }

  dataframe = temp_df_to_build

  return(dataframe)

}
