package ara

import java.io.File
import java.io.PrintWriter
import scala.io.Source
import ara.Cluster._
import edu.northwestern.at.utils.math.statistics.Descriptive

object AraUtilities extends MTBCclusters {

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
      def medianCov(c: String): Double = {
        val medianCovMap = markers.map(_ match {
          case (lineage, linMarkers) => {
            val arr = linMarkers.filter(_.isPresent).map(_.count.toDouble).sorted.toArray
            (lineage, Descriptive.median(arr))
          }
        })
        if (c.hasZeroMarkers || c == "MTBC") 0 else medianCovMap(c)
      }

      /** Total median coverage at MTBC root */
      val totalMedianCov = medianCov("L1-L2-L3-L4") + medianCov("L5-L6-LB")

      /**
       * Presence of cluster based on number detected markers
       *  Clusters with 0 markers could be present, so function return true for these clusters.
       */
      def presence(c: String): Boolean = {
        if (c.hasZeroMarkers || c == "MTBC") return true
        val totalMarkers = totalMarkersPerLineage(c)
        val presentMarkers = linCountsPresent(c)
        if (c.hasReference) {
          presentMarkers > (totalMarkers / 30)
        } else {
          presentMarkers > (totalMarkers / 10)
        }
      }

      /** Recursive call to get tree path from leaf to root */
      def getPath(str: String): List[String] = {
        def recursive(ls: List[String]): List[String] = ls.head match {
          case c if (c.hasAncestor) => recursive(c.getAncestor :: ls)
          case _ => ls
        }
        recursive(List(str))
      }


      /*def intersectAll(ls: List[List[String]]): List[String] = {
        if (ls.isEmpty) return List.empty[String]
        def recursiveIntersect(ls1: List[String], ls2: List[List[String]]): List[String] = ls2 match {
          case head :: tail => recursiveIntersect(ls1 intersect head, tail)
          case Nil => ls1
        }
        recursiveIntersect(ls.head, ls.tail)
      }*/

      def medianCovPath(ls: List[String]): Double = {
        val covPerCluster = ls.filterNot(_.hasZeroMarkers).map(c => medianCov(c)).sorted.toArray
        Descriptive.median(covPerCluster)
      }

      def meanCovPath(ls: List[String]): Double = {
        val covPerCluster = ls.filterNot(_.hasZeroMarkers).map(c => medianCov(c))
        if (covPerCluster.isEmpty) 0
        else (covPerCluster).foldLeft(0.toDouble)(_ + _) / covPerCluster.size
      }

      /**
       *  Interpret results 
       */
      val presentLeaves = clusters.filter(_.isLeafCluster).filter { c => presence(c) }
      presentLeaves.sorted.foreach(c => if (c.hasZeroMarkers) println(c + "\t" + 0 + "\t" + 0) else println(c + "\t" + averageCov(c) + "\t" + medianCov(c)))
      println

      val paths = presentLeaves.map(getPath(_)).filterNot { p =>
        val presencePath = p.map(c => presence(c))
        presencePath.contains(false)
      } // Filter out paths that have gaps, clusters in path for which half of the markers is not detected.

      /*.map{path =>
          val cov = path.filterNot(_.hasZeroMarkers).map(medianCov(_))
          (path zip cov)
      }*/
      println("Present leave(s)")
      presentLeaves.foreach(println)
      println
      println("Present path(s)")
      paths.foreach(println)
      println
      
      /** 
       *  Print output
       */
      val pw = new PrintWriter(config.output)

      pw.println("# Ara Results")
      pw.println("# Command: " + args.mkString(" "))
      pw.println("# ")
      pw.println("Predicted groups: " + paths.map(_.last).mkString(", "))
      pw.println("Mixed infection: " + (paths.size > 1))

      
      val strainCountPerCluster = paths.flatten.groupBy(identity).mapValues(_.size)
      val presentStrains = strainCountPerCluster.map(_._2).max

      pw.println(presentStrains + " present strains/paths")

      val separatedPaths = (1 to presentStrains).toList.map { number =>
        val path = paths.map(_.filter(strainCountPerCluster(_) == number)).distinct.filterNot(_.isEmpty)
        (number -> path)
      }.toMap
      separatedPaths.foreach(println)
      println

      val meanCovPerPath = separatedPaths.map {
        _ match {
          case (number, pathList) => (number, pathList.map(path => (meanCovPath(path), medianCovPath(path), path)))
        }
      }
      meanCovPerPath.foreach(println)
      println

      val rootNodes = meanCovPerPath.flatMap(_._2).filter(p => p._3.head == "L1-L2-L3-L4" || p._3.head == "L5-L6-LB" || p._3.head == "MTBC")
      val rootMeanCov = rootNodes.map(_._1).foldLeft(0.toDouble)(_ + _)
      val rootMedianCov = rootNodes.map(_._2).foldLeft(0.toDouble)(_ + _)
      rootNodes.foreach(println)
      println("Mean coverage root node: " + rootMeanCov + " (1)")
      println("Median coverage root node: " + rootMedianCov + " (1)")

      if (meanCovPerPath.size > 1) { // Mixed infection, estimate frequencies

        val covL1toL4 = meanCovPerPath.flatMap(_._2).filter(p => p._3.head == "L1-L2-L3-L4")
        println("Path coverage L1-L2-L3-L4: " + covL1toL4)
        val covL5toLB = meanCovPerPath.flatMap(_._2).filter(p => p._3.head == "L5-L6-LB")
        println("Path coverage LL5-L6-LB: " + covL5toLB)
        println
        
        
        
      } else { // Not a mixed infection
        println("One present strain/path")
        println(paths.flatten.map(c => (c, medianCov(c))))
      }

      pw.close

    }

  }
}