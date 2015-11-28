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

      /*def medianCovPath(ls: List[String]): Double = {
        val covPerCluster = ls.filterNot(_.hasZeroMarkers).map(c => medianCov(c)).sorted.toArray
        Descriptive.median(covPerCluster)
      }*/

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
      paths.foreach(p => println(p.map(c => (medianCov(c), c))))
      println

      /**
       *  Print output
       */
      val pw = new PrintWriter(config.output)

      pw.println("# Ara Results")
      pw.println("# Command: " + args.mkString(" "))
      pw.println("# ")

      val pathCount = paths.flatten.groupBy(identity).mapValues(_.size)
      val presentStrains = pathCount.map(_._2).max
      pw.println("# " + presentStrains + " present strain(s)/path(s)\n")

      val separatedPaths = (1 to presentStrains).toList.map { number =>
        val path = paths.map(_.filter(pathCount(_) == number)).distinct.filterNot(_.isEmpty)
        (number -> path)
      }.toMap
      println("Separate path(s)")
      separatedPaths.foreach(println)
      println

      val splitPaths = separatedPaths.flatMap(_._2).toList

      val meanCovPerPath = splitPaths.map { path =>
        (meanCovPath(path), path)
      }
      println("Mean coverage per path")
      meanCovPerPath.foreach(println)
      println

      if (presentStrains > 1) { // Mixed infection, estimate frequencies

        def getLevel(p: List[String]): Int = p match {
          case head :: tail => head match {
            case "MTBC" => 0
            case x => { 1 + getLevel(splitPaths.filter(_.contains(x.getAncestor)).head) }
          }
          case Nil => 0
        }

        val pathNumbers = meanCovPerPath.map {
          _ match {
            case (depth, path) => {
              if (path.head == "MTBC") {
                if (path.size == 1) {
                  val rootNodes = meanCovPerPath.filter(p => p._2.head == "L1-L2-L3-L4" || p._2.head == "L5-L6-LB" || p._2.head == "MTBC")
                  val rootMeanCov = rootNodes.map(_._1).foldLeft(0.toDouble)(_ + _)
                  (1.toDouble, getLevel(path), rootMeanCov, rootMeanCov, path)
                } else (1.toDouble, getLevel(path), depth, depth, path)
              } else {
                val sibling = path.head.getSibling
                val siblingPath = meanCovPerPath.filter(_._2.contains(sibling)).head
                val totalCovSiblings = siblingPath._1 + depth
                val frequency = depth / totalCovSiblings
                (frequency, getLevel(path), totalCovSiblings, depth, path)
              }
            }
          }
        }.sortBy(_._2)
        println("Frequency, level, totalDepth, depth, path")
        pathNumbers.foreach(println)
        println

        def getAncestralNodes(p: List[String]): List[String] = p match {
          case head :: tail => {
            if (head == "MTBC") List(head)
            else head :: getAncestralNodes(splitPaths.filter(_.contains(head.getAncestor)).head)
          }
          case Nil => p
        }
        
        val frequencies = pathNumbers.filter(p => p._5.map(_.isLeafCluster).contains(true)).map{p =>
          val ancestors = getAncestralNodes(p._5)
          println(ancestors)
          val freqArr = ancestors.map(a => pathNumbers.filter(_._5.contains(a)).head).map(_._1)
          println(freqArr)
          println
          (p._5.last, (freqArr.foldLeft(1.toDouble)(_ * _)))
        }
        println("Frequencies")
        frequencies.foreach(println)
        println
        
        /**
         * Print Output
         */
        
        pw.println("Predicted group(s): " + frequencies.map(p => p._1 + "(" + p._2 + ")").mkString(", "))
        def printNumbers(ls: List[(Double, Int, Double, Double, List[String])]) = {
          pw.println("Mixed infection: TRUE\n")
          ls.foreach {
            _ match {
              case (freq, lvl, siblingsDepth, depth, path) => {
                val space = "\t" * lvl
                pw.println(space + path.mkString(" -> "))
                pw.println(space + "Mean read depth: " + depth)
                pw.println(space + "Frequency estimate: " + freq + " (" + depth + "/" + siblingsDepth + ")")
                pw.println
              }
            }
          }
        }
        printNumbers(pathNumbers)

      } else if (presentStrains == 1) { // Not a mixed infection

        pw.println("Predicted group(s): " + paths.map(_.last).mkString(", "))
        pw.println("Mixed infection: FALSE\n")
        pw.println(paths.head.mkString(" -> "))
        pw.println("Mean read depth: " + meanCovPerPath.head._1)

      } else { // No path detected

        pw.println("Predicted groups: NONE")
        pw.println("Mixed infection: FALSE")

      }

      pw.close

    }

  }
}