package ara

import java.io.File
import java.io.PrintWriter
import scala.io.Source
import ara.Cluster._
import edu.northwestern.at.utils.math.statistics.Descriptive
import java.util.Calendar

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

      /** Count total number of markers per cluster. */
      val totalMarkersPerLineage = cNumbers.mapValues(_("markers")) //markers.mapValues(_.size)
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

      /**
       * Presence of cluster based on number detected markers
       * Clusters with 0 markers could be present, so function return true for these clusters.
       */
      def presence(c: String): Boolean = {
        if (!c.isCluster) false
        else if (c.hasZeroMarkers || c == "MTBC") true
        else {
          val totalMarkers = totalMarkersPerLineage(c)
          val presentMarkers = linCountsPresent(c)
          if (c.hasReference) {
            presentMarkers > (totalMarkers / 24)
          } else {
            presentMarkers > (totalMarkers / 8)
          }
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

      /** Mean coverage of of clusters in path */
      def meanCovPath(ls: List[String]): Double = {
        val covPerCluster = ls.filterNot(_.hasZeroMarkers).map(c => medianCov(c))
        if (covPerCluster.isEmpty) 0
        else (covPerCluster).foldLeft(0.toDouble)(_ + _) / covPerCluster.size
      }

      /** Split paths by present strains */
      def splitPaths(ls: List[List[String]]): List[List[String]] = {
        if (ls.isEmpty) List.empty[List[String]]
        else {
          val pathCount = ls.flatten.groupBy(identity).mapValues(_.size)
          val maxPaths = pathCount.map(_._2).max
          val separatedPaths = (1 to maxPaths).toList.map { number =>
            val path = ls.map(_.filter(pathCount(_) == number)).distinct.filterNot(_.isEmpty)
            (number -> path)
          }.toMap
          separatedPaths.flatMap(_._2).toList
        }

      }

      def meanCovSeparatePaths(ls: List[List[String]]): List[(Double, List[String])] = {
        val res = splitPaths(ls).map { path =>
          (meanCovPath(path), path)
        }
        res
      }

      def childIsPresent(c: String): Boolean = {
        val children = c.children
        if (children(0).isCluster && children(1).isCluster) (presence(children(0)) || presence(children(1)))
        else false
      }

      def ancestorIsPresent(c: String): Boolean = {
        val ancestor = c.getAncestor
        presence(c)
      }

      /** Check possible presence of (separate) paths where each cluster has 0 markers */
      def pathWith0Markers(p: List[String]): Boolean = {
        !p.map(c => c.hasZeroMarkers).contains(false)
      }

      def printCovSeparatePaths(ls : List[(Double, List[String])]) : Unit = {
        println("Mean coverage per separate path (markers, cluster)")
        ls.foreach(p => println(p._1 + "\t" + p._2.map { c =>
          val markers = if (linCountsPresent.contains(c)) linCountsPresent(c) + "/" + (if (c.hasReference) totalMarkersPerLineage(c) / 3 else totalMarkersPerLineage(c)) else "0"
          (markers, c)
        }))
        println
      }
      
      /** Remove absent paths */
      def removeAbsent(ls: List[List[String]]): List[List[String]] = {
        val meanCovPerPath = meanCovSeparatePaths(ls)
        printCovSeparatePaths(meanCovPerPath)
        def listAbsent(ls: List[(Double, List[String])]): List[String] = {
          val absentClusters = ls.filter {
            _ match {
              case (mean, path) => {

                if (path.equals(List("MTBC"))) {
                  false
                } else {
                  val sibling = path.head.getSibling
                  if (sibling != null) { // sibling path exists
                    val siblingPath = ls.filter(_._2.contains(sibling)).head
                    val siblingCov = siblingPath._1
                    val ancestor = path.head.getAncestor
                    val ancestorCov = ls.filter(_._2.contains(ancestor)).head._1
                    if (pathWith0Markers(path)) { // sibling with no markers
                      if (((ancestorCov - siblingCov) > 25)) { // If there is missing read depth >25 compared to ancestor, then path with 0 markers is present
                        false
                      } else { // No missing read depth compared to ancestor, so path with 0 markers is absent
                        true
                      }
                    } else { // sibling with markers
                      if (ancestor == "MTBC") {
                        false
                      } else {
                        val totalCovSiblings = siblingCov + mean
                        if ((totalCovSiblings / ancestorCov) > 1.8) { // If unbalanced cluster sizes of siblings
                          println("path: " + path + ", sibling: " + sibling)
                          println(totalCovSiblings + " / " + ancestorCov + " = " + (totalCovSiblings / ancestorCov))
                          if (path.head.getSize > sibling.getSize) {
                            println(path.head.getSize + " > " + sibling.getSize)
                            true
                          } else {
                            false
                          }
                        } else {
                          false
                        }
                      }

                    }
                  } else { // sibling path not exists                  
                    false
                  }
                }
              }
            }
          }.flatMap(_._2)
          absentClusters
        }
        val absent = listAbsent(meanCovPerPath)

        if (absent.isEmpty) {
          ls
        } else {
          println("Absent cluster with 0 markers or unbalanced cluster")
          absent.foreach(println)
          println
          val presentPaths = ls.filter { p =>
            p.intersect(absent).isEmpty
          }
          removeAbsent(presentPaths)
        }
      }

      /** Get numbers for separate paths */
      def pathNums(ls: List[(Double, List[String])]): List[(Double, Int, Double, Double, List[String])] = {
        def getLevel(p: List[String]): Int = p match {
          case head :: tail => head match {
            case "MTBC" => 0
            case x => 1 + getLevel(ls.map(_._2).filter(_.contains(x.getAncestor)).head)
          }
          case Nil => 0
        }

        val res = ls.map {
          _ match {
            case (depth, path) => {
              if (path.head == "MTBC") { // root path
                if (path.size == 1) { // root path only contains cluster MTBC
                  val rootNodes = ls.filter(p => p._2.head == "L1-L2-L3-L4" || p._2.head == "L5-L6-LB" || p._2.head == "MTBC")
                  val rootMeanCov = rootNodes.map(_._1).foldLeft(0.toDouble)(_ + _)
                  (1.toDouble, getLevel(path), rootMeanCov, rootMeanCov, path)
                } else (1.toDouble, getLevel(path), depth, depth, path)
              } else { //non-root path
                val sibling = path.head.getSibling
                val siblingPath = ls.filter(_._2.contains(sibling)).head
                val ancestorCov = ls.filter(_._2.contains(path.head.getAncestor)).head._1
                if (pathWith0Markers(siblingPath._2)) { // Misses coverage compared to ancestor (> 25 reads)
                  val frequency = depth / ancestorCov
                  (frequency, getLevel(path), ancestorCov, depth, path)
                } else if (pathWith0Markers(path)) { // Present path with 0 markers
                  val frequency = 1 - (siblingPath._1 / ancestorCov)
                  (frequency, getLevel(path), ancestorCov, depth, path)
                } else {
                  val totalCovSiblings = siblingPath._1 + depth
                  val frequency = depth / totalCovSiblings
                  (frequency, getLevel(path), totalCovSiblings, depth, path)
                }
              }
            }
          }
        }.sortBy(_._2)
        res
      }

      /** Multiply frequencies of separate paths to ancestor */
      def getFrequencies(ls: List[(Double, Int, Double, Double, List[String])]): List[(String, Double)] = {
        /** Get ancestral paths for separate path */
        def getAncestralNodes(p: List[String]): List[String] = p match {
          case head :: tail => head match {
            case "MTBC" => List(head)
            case x => x :: getAncestralNodes(ls.map(_._5).filter(_.contains(x.getAncestor)).head)
          }
          case Nil => Nil
        }

        val res = ls.filterNot(p => (p._5.last.children.map(cc => presence(cc)).contains(true))).map { p =>
          val ancestors = getAncestralNodes(p._5)
          val freqArr = ancestors.map(a => ls.filter(_._5.contains(a)).head).map(_._1)
          (p._5.last, (freqArr.foldLeft(1.toDouble)(_ * _)))
        }
        res
      }
      
      def removeLowFreqPaths(ls: List[List[String]]): List[List[String]] = {
        val meanCovPerPath = meanCovSeparatePaths(ls)
        printCovSeparatePaths(meanCovPerPath)
        val pathNumbers = pathNums(meanCovPerPath)
        val frequencies = getFrequencies(pathNumbers)
        val lowFreqStrains = frequencies.filter(_._2 < 0.01)

        if (lowFreqStrains.isEmpty) {
          ls
        } else {
          println("Low frequency strain(s)/path(s)")
          lowFreqStrains.foreach(println)
          println
          val present = ls.filter(p => p.intersect(lowFreqStrains.map(_._1)).isEmpty)
          removeLowFreqPaths(present)
        }
      }

      /**
       *  Interpret results and print results in console
       */

      val clusters = mtbcClusters.filterNot(_ == "MTBC").filter { c => presence(c) }.filterNot(c => childIsPresent(c)).filter(ancestorIsPresent(_))
      println("Possible present end cluster(s)\tMean coverage\tMedian coverage")
      clusters.sorted.foreach(c => if (c.hasZeroMarkers) println(c + "\t" + 0 + "\t" + 0) else println(c + "\t" + averageCov(c) + "\t" + medianCov(c)))
      println

      val paths = clusters.map(getPath(_)).filterNot { p =>
        val path = p.map(c => presence(c))
        path.contains(false)
      }.filterNot(_.last == "L1-L2-L3-L4").filterNot(_.last == "L5-L6-LB")
      println("Possible present path(s)")
      paths.foreach(p => println(p.map(c => (medianCov(c), c))))
      println

      val filteredPaths = removeAbsent(paths)
      println("Filtered path(s)")
      filteredPaths.foreach(println)
      println
      
      val filteredPaths2 = removeLowFreqPaths(filteredPaths)
      println("Filtered path(s)")
      filteredPaths2.foreach(println)
      println
      
      val separatePathsCov = meanCovSeparatePaths(filteredPaths2)
      
      val pathNumbers = pathNums(separatePathsCov)
      println("Frequency, level, totalDepth, depth, path")
      pathNumbers.foreach(println)
      println

      /**
       *  Print output
       */
      val pw = new PrintWriter(config.output)
      pw.println("# Ara Results")
      pw.println("# Command: " + args.mkString(" "))
      pw.println("# Date: " + Calendar.getInstance.getTime)
      val strains = pathNumbers.filterNot(p => childIsPresent(p._5.last)).size
      pw.println("# " + strains + " present strain(s)/path(s)")
      pw.println

      def printNumbers(ls: List[(Double, Int, Double, Double, List[String])]) = {
        if (ls.size < 2) pw.println("Mixed sample: FALSE\n")
        else pw.println("Mixed sample: TRUE\n")
        ls.foreach {
          _ match {
            case (freq, lvl, siblingsDepth, depth, path) => {
              val space = "\t" * lvl
              pw.println(space + (if (path.size == 1) path.head else (path.head + " -> " + path.last)))
              pw.println(space + "Mean read depth: " + depth)
              pw.println(space + "Frequency estimate: " + freq + " (" + depth + "/" + siblingsDepth + ")")
              pw.println
            }
          }
        }
      }

      if (pathNumbers.size > 1) { // Mixed infection, estimate frequencies
        val frequencies = getFrequencies(pathNumbers)
        println("Total frequencies")
        frequencies.foreach(println)
        println
        pw.println("Predicted group(s): " + frequencies.map(p => p._1 + "(" + p._2 + ")").mkString(", "))
      } else if (pathNumbers.size == 1) { // Not a mixed infection
        val path = pathNumbers.head._5
        val cov = pathNumbers.head._4
        pw.println("Predicted group(s): " + path.last)
      } else { // No path detected
        pw.println("Predicted group(s): NONE")
      }

      printNumbers(pathNumbers)

      pw.close

    }

  }
}