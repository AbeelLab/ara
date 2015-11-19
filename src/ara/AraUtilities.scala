package ara

import java.io.File
import java.io.PrintWriter
import scala.io.Source
import ara.Cluster._

object AraUtilities {

  case class Config(
    val markers: File = null,
    val output: String = null)

  def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar ara.jar interpret") {
      opt[File]('m', "markers") required () action { (x, c) => c.copy(markers = x) } text ("Input file with mapped markers.")
      opt[String]('o', "output") required () action { (x, c) => c.copy(output = x + ".interpret.ara") } text ("Output name.")
    }

    parser.parse(args, Config()) map { config =>

      val pw = new PrintWriter(config.output)

      val markers = Source.fromFile(config.markers).getLines.filterNot(_.startsWith("#")).toList.dropRight(2).map {
        _ match {
          case Marker(m, c, p) => new ClusterMarker(m, c, p)
        }
      }.groupBy(_.lineage)
      //markers.foreach(pw.println)
      //println(markers.size  + " clusters")
      val clusters = markers.keysIterator.toList
      //clusters.foreach(println)

      /** Count total number of markers per cluster. */
      val totalMarkersPerLineage = markers.mapValues(_.size)
      //totalMarkersPerLineage.toList.sortBy(_._1).foreach(println)
      //println("Total markers per lineage")
      //println(totalMarkersPerLineage.foldLeft(0)(_ + _._2))

      /** Total number of mapped reads per cluster. */
      val markerCounts = markers.map(_ match {
        case (lineage, linMarkers) => (lineage, linMarkers.map(_.count).foldLeft(0)(_ + _))
      })
      //println("markerCounts")
      //markerCounts.toList.sortBy(_._1).foreach(println)

      /** Number of present markers per cluster. */
      val linCountsPresent = markers.map(_ match {
        case (lineage, linMarkers) => (lineage, linMarkers.filter(_.isPresent).size)
      })
      //println("LinCountsPresent")
      //linCountsPresent.foreach(pw.println)

      /** Average read depth of present markers per cluster */
      val averageCov = markers.filter(_ match {
        case (lineage, linMarkers) => linCountsPresent(lineage) > 0
      }).map(_ match {
        case (lineage, linMarkers) => (lineage, linMarkers.foldLeft(0)(_ + _.count) / linCountsPresent(lineage))
      })
      //averageCov.foreach(println)

      pw.println("$$\tTotal mapped reads\tTotal markers\tPresent markers\tMean read-depth/present marker")

      clusters.sorted.foreach { c =>
        val avgCov = if (averageCov.contains(c)) averageCov(c) else 0
        pw.println(c + "\t" + markerCounts(c) + "\t" + totalMarkersPerLineage(c) + "\t" + linCountsPresent(c) + "\t" + avgCov)
      }

      pw.close

      def getPresentClusters(): Unit = {
        /**val rootClusters = List("L1-L2-L3-L4", "L5-L6-LB").map { c =>
          val snpPositions = if (c.hasReference) totalMarkersPerLineage(c) / 3 else totalMarkersPerLineage(c)
          val present = linCountsPresent(c)
          if (present > (snpPositions / 2)) (c, true) else (c, false)
        }*/
        
        /** Return presence is true if more than 50% of SNP positions has been detected. */
        def getPresence(c: String): Boolean = {
          val snpPositions = if (c.hasReference) totalMarkersPerLineage(c) / 3 else totalMarkersPerLineage(c)
          val present = linCountsPresent(c)
          if (present > (snpPositions / 2)) true else false
        }
        
        val path1 = List("L1-L2-L3-L4", "L1")
        val path2 = List("L1-L2-L3-L4", "L2-L3-L4", "L2-L3", "L2")
        val path3 = List("L1-L2-L3-L4", "L2-L3-L4", "L2-L3", "L3")
        val path4 = List("L1-L2-L3-L4", "L2-L3-L4", "L4")
        val path5 = List("L5-L6-LB", "L5")
        val path6 = List("L5-L6-LB", "L6-LB", "L6")
        val pathB = List("L5-L6-LB", "L6-LB", "LB")
        
        
        
        /**rootClusters.foreach(_._2 match {
          case true => 
          case false =>
        })
        rootClusters.foreach(println) */
      }

      getPresentClusters()

    }

  }
}