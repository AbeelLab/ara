package ara

import java.io.File
import java.io.PrintWriter
import scala.io.Source
import ara.Cluster._
import edu.northwestern.at.utils.math.statistics.Descriptive

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

      /** Median read depth of present markers per cluster */
      val medianCov = markers.map(_ match {
        case (lineage, linMarkers) => {
          val arr = linMarkers.filter(_.isPresent).map(_.count.toDouble).sorted.toArray
          (lineage, Descriptive.median(arr))
        }
      })
      
      /** Total median coverage at MTBC root */
      val totalMedianCov = medianCov("L1-L2-L3-L4") + medianCov("L5-L6-LB")

      
      val presentClusters = clusters.filter { c =>
        val totalMarkers = totalMarkersPerLineage(c)
        val presentMarkers = linCountsPresent(c)
        if (c.hasReference) {
          presentMarkers > (totalMarkers / 6)
        } else {
          presentMarkers > (totalMarkers / 2)
        }
      }
      //val filtered = presentClusters.filter(c => presentClusters.contains(c.getAncestor) || c.getAncestor == "MTBC")
      //
      presentClusters.sorted.foreach(c => println(c + "\t" + averageCov(c) + "\t" + medianCov(c)))
      println
    

      /** Recursive call to get tree path from leaf to root */
      def getPath(str: String): List[String] = {
        def recursive(ls: List[String]): List[String] = ls.head match {
          case c if (c.hasAncestor) => recursive(c.getAncestor :: ls)
          case _ => ls
        }
        recursive(List(str))
      }
      
      val leaves = presentClusters.filter(_.isLeafCluster).map(getPath(_)).map{path =>
          val cov = path.filterNot(_.hasZeroMarkers).map(medianCov(_))
          (path zip cov)
      }
      println("Leaves")
      leaves.foreach(l => println(l))
      println
      println("paths")
      val paths = leaves.filterNot(c => c.map(_._2).contains(0))//Filter out leaf cluster if it has an ancestor with 0 read depth
      paths.foreach(println)
      println
      val intersect = paths.flatten.distinct
      //intersect.foreach(println)
      val overlappingPath = paths.flatten.groupBy(identity).mapValues(_.size).filter(c => c._2 > 1).keysIterator.toList.sorted
      println("Overlapping path")
      println(overlappingPath)
      val splitPaths = paths.map(_.filterNot(c => overlappingPath.map(_._1).contains(c._1)))
      println("Split paths")
      splitPaths.foreach(println)
      
      
      /** Print output*/
      val pw = new PrintWriter(config.output)      
      /**pw.println("$$\tTotal mapped reads\tTotal markers\tPresent markers\tMean read-depth/present marker")
      clusters.sorted.foreach { c =>
        val avgCov = if (averageCov.contains(c)) averageCov(c) else 0
        pw.println(c + "\t" + markerCounts(c) + "\t" + totalMarkersPerLineage(c) + "\t" + linCountsPresent(c) + "\t" + avgCov)
      }**/
      pw.println("# Ara Results")
      pw.println("# Command: " + args.mkString(" "))
      pw.println("# ")
      pw.println("Predicted groups: " + paths.map(_.last).mkString(", "))
      pw.println("Mixed infection: " + (paths.size > 1))
      pw.println
      intersect.foreach(pw.println)
      
      pw.close

      
      

    }

  }
}