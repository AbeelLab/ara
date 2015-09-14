package macaw2

import java.io.File
import scala.io.Source
import java.io.PrintWriter
import edu.northwestern.at.utils.math.statistics.FishersExactTest


object MacawUtilities {
  
  case class Config(
    val snpTyperOutput: String = null, 
    val result: String = null
  )
  
  
  def main(args: Array[String]) {
    
    val parser = new scopt.OptionParser[Config]("java -jar ara.jar interpret-GT") {
      opt[String]('m', "markers") required() action { (x, c) => c.copy(snpTyperOutput = x) } text ("Output file of MacawSNPTyper.")
      opt[String]('o', "output") required() action { (x, c) => c.copy(result = x + ".interpret-MI.ara") } text ("Output name for the file with results.")
    }
    
    parser.parse(args, Config()) map { config =>
      
      /** Read markers from file and group by cluster. */
      val markers = Source.fromFile(new File(config.snpTyperOutput)).getLines.filterNot(_.startsWith("#")).toList.dropRight(2).map { _ match {
        case Marker(m, c, p) => new ClusterMarker(m, c, p)
      }}.groupBy(_.lineage)
      
      /** Count total number of markers per cluster. */
      val totalMarkersPerLineage = markers.mapValues(_.size)
      totalMarkersPerLineage.foreach(println)
      println(totalMarkersPerLineage.foldLeft(0)(_ + _._2))
      
      /** Total number of mapped reads per cluster. */
      val markerCounts = markers.map(_ match {
        case (lineage, linMarkers) => (lineage, linMarkers.map(_.count).foldLeft(0)(_ + _))
      })
      println("markerCounts")
      markerCounts.foreach(println)
      
      /** Number of present markers per cluster. */
      val linCountsPresent = markers.map(_ match {
        case (lineage, linMarkers) => (lineage, linMarkers.filter(_.isPresent).size)
      })
      println("LinCountsPresent")
      linCountsPresent.foreach(println)
      
      /** Number of absent markers per cluster. */
      val linCountsAbsent = markers.map(_ match {
        case (lineage, linMarkers) => (lineage, linMarkers.filterNot(_.isPresent).size)
      })
      println("LinCountsAbsent")
      linCountsAbsent.foreach(println)
      
      /** Sum number of present and absent markers over all clusters. */
      val sumAllPresent = linCountsPresent.foldLeft(0)(_ + _._2)
      val sumAllAbsent = linCountsAbsent.foldLeft(0)(_ + _._2)
      
      println("sumAllPresent: " + sumAllPresent)
      println("sumAllAbsent: " + sumAllAbsent)
      
      /** List of clusters to detect. */
      val lineages = markers.keysIterator.toList.sorted
        
      /** Calculate fisher p-values for each cluster. */
      val fisherPValues = lineages.map{lin =>
          val pValue = FishersExactTest.fishersExactTest(linCountsPresent(lin), sumAllPresent - linCountsPresent(lin), linCountsAbsent(lin), sumAllAbsent - linCountsAbsent(lin))(2) * lineages.size
          if (pValue > 1) (lin, 1.toFloat)
          else (lin, pValue.toFloat)        
      }.toMap
      println("fisherPValues: ")
      fisherPValues.foreach(println)
      
      /** Average read depth per cluster: # of mapped reads per clusters / # of markers per cluster. */
      val avgMarkerPerLin = lineages.map(lin => (lin, markerCounts(lin).toFloat / totalMarkersPerLineage(lin))).toMap
      println("avgMarkerPerLin: ")
      avgMarkerPerLin.foreach(println)
      
      /** Percentage of present clusters. */
      val allLinDepth = avgMarkerPerLin.foldLeft(0.toFloat)(_ + _._2)
      val linPercentTotal = lineages.map(lin => (lin, avgMarkerPerLin(lin) / allLinDepth)).toMap
      
      
      /** Detected clusters */
      val predictedLin = fisherPValues.filter( _ match {
        case (lin, pValue) => pValue < 0.05
      }).map(_._1)
      print("predictedLin: " + predictedLin.mkString(", "))
      
      
      /**
       * Print results to file
       */
      val pw = new PrintWriter(new File(config.result))
      pw.println("# Results")      
      predictedLin.size match {
        case x if (x > 1) => {
          val totalLinDepth = predictedLin.foldLeft(0.toFloat)(_ + avgMarkerPerLin(_))
          val linPercentEnd = predictedLin.map(lin => lin + " (" + (avgMarkerPerLin(lin) / totalLinDepth) + ")").mkString(", ")
          pw.println("Predicted group(s):\t" + linPercentEnd + "\n" + "Mixed sample?:\tTRUE\n")
        }
        case 0 => pw.println("Predicted group(s):\tNONE" + "\n" + "Mixed sample?:\tFALSE\n")
        case 1 => pw.println("Predicted group(s):\t" + predictedLin.head + "\n" + "Mixed sample?:\tFALSE\n")        
      }
      pw.println("# Details")
      pw.println(lineages.map{lin => 
        (lin + " Present markers:\t" + linCountsPresent(lin) + "\n" + 
            lin + " Absent markers:\t" + linCountsAbsent(lin) + "\n" + 
            lin + " Raw marker total:\t" + markerCounts(lin) + " (" + linPercentTotal(lin) + ")" + "\n" + 
            lin + " Corrected p-value:\t" + fisherPValues(lin))
      }.mkString("\n\n"))
      
      pw.close

    }    
  }
}