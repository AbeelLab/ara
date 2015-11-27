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
      } 
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
      

      
      val strainCountPerCluster = paths.flatten.groupBy(identity).mapValues(_.size)
      val presentStrains = strainCountPerCluster.map(_._2).max
      pw.println("# " + presentStrains + " present strains/paths\n")

      val separatedPaths = (1 to presentStrains).toList.map { number =>
        val path = paths.map(_.filter(strainCountPerCluster(_) == number)).distinct.filterNot(_.isEmpty)
        (number -> path)
      }.toMap
      separatedPaths.foreach(println)
      println

      val meanCovPerPath = separatedPaths.map {
        _ match {
          case (number, pathList) => (number, pathList.map(path => (meanCovPath(path), path)))
        }
      }
      meanCovPerPath.foreach(println)
      println

      val rootNodes = meanCovPerPath.flatMap(_._2).filter(p => p._2.head == "L1-L2-L3-L4" || p._2.head == "L5-L6-LB" || p._2.head == "MTBC")
      val rootMeanCov = rootNodes.map(_._1).foldLeft(0.toDouble)(_ + _)
      //val rootMedianCov = rootNodes.map(_._2).foldLeft(0.toDouble)(_ + _)
      rootNodes.foreach(println)
      println("Mean coverage root node: " + rootMeanCov + " (1)")
      //println("Median coverage root node: " + rootMedianCov + " (1)")

      
      if (presentStrains > 1) { // Mixed infection, estimate frequencies
        
        val freqs = meanCovPerPath(1).map(p => ((p._1 / rootMeanCov), p._1, p._2.last))
        pw.println("Predicted groups: " + freqs.map(p => p._3 + "(" + p._1 +")").mkString(", "))
        pw.println("Mixed infection: TRUE\n")
        
        val rootPath = meanCovPerPath(presentStrains).head // Path with all strains
        pw.println("Root path: " + rootPath._2.mkString(" -> "))
        pw.println("Mean read coverage: " + rootMeanCov)
        pw.println("Frequency estimate: " + 1.0)
        pw.println
        

        for (number <- (1 to presentStrains - 1).reverse){
          println(number)
          val path = meanCovPerPath(number).map{p => 
            val ancestor = p._2.head.getAncestor
            println(ancestor)
            val ancestorPath = meanCovPerPath.flatMap(_._2).filter(_._2.contains(ancestor)).head
            
            val space = "\t" * (presentStrains - number)
            pw.println(space + p._2.mkString(" -> "))
            pw.println(space + "Mean read coverage: " + p._1)
            if (ancestor == "MTBC") pw.println(space + "Frequency estimate: " + (p._1 / rootMeanCov))
            else pw.println(space + "Frequency estimate: " + (p._1 / ancestorPath._1))
            pw.println
          }
          
        }
                
      } else if (meanCovPerPath.size == 1){ // Not a mixed infection
        
        pw.println("Predicted groups: " + paths.map(_.last).mkString(", "))
        pw.println("Mixed infection: FALSE")
        
        println("One present strain/path")
        println(paths.flatten.map(c => (c, medianCov(c))))
        
        
      } else { // No path detected
        pw.println("Predicted groups: NONE")
        pw.println("Mixed infection: FALSE")
      }

      pw.close

    }

  }
}