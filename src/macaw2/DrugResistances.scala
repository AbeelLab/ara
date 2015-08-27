package macaw2

import java.io.File
import scala.io.Source
import java.io.PrintWriter


object DrugResistances {
  
  case class Config(
    val snpTyperOutput: String = null, 
    val result: String = null
  )
  
  class Marker(val mutInfo: String, val count: Int, val isPresent: Boolean){
    override def toString(): String = mutInfo + "\t" + count + "\t" + (if (isPresent) "present" else "absent")
    val arr = mutInfo.split("_")
    val coordinate = arr(0).drop(1).dropRight(1).toInt
    val markerType = arr(1)
    val drugs = arr(2)
  }

  object Marker {
    def unapply(s: String): Option[(String, Int, Boolean)] = {
      val line = s.split("\t")
      val mutInfo = line(0)
      val count = line(1).toInt
      val presence = if (line(3) == "present") true else false
      Some(mutInfo, count, presence)
    }
  }
  
  def main(args: Array[String]){
    
    val parser = new scopt.OptionParser[Config]("java -jar Ara.jar interpret-DR"){
      opt[String]('m', "markers") required() action { (x, c) => c.copy(snpTyperOutput = x) } text ("Output file of MacawSNPTyper.")
      opt[String]('o', "output") required() action { (x, c) => c.copy(result = x) } text ("Output name for the file with results.")
    }
    
    class findInList(val ls: List[Marker]){
      val presentMarkers = ls.filter(_.isPresent)
      def hasResistance = {
        if (presentMarkers.exists(_.markerType == "resistance")) true
        else false
      }
      def hasSusceptibility = {
        if (presentMarkers.exists(_.markerType == "susceptibility")) true
        else false
      }
    }
    
    implicit def checkMarkers(ls: List[Marker]) = new findInList(ls)
    
    parser.parse(args, Config()) map { config =>
      val snpTypes = new File(config.snpTyperOutput)
      val outFile = new File(config.result)
      
      val markers = Source.fromFile(snpTypes).getLines.filterNot(_.startsWith("#")).toList.dropRight(2).map {_ match {
          case Marker(mutInfo, count, isPresent) => new Marker(mutInfo, count, isPresent)
      }}.groupBy(_.coordinate)//.map(_._2).toList.sortBy(_(1).coordinate)
      
      //markers.foreach(println)
      println(markers.size + " genome positions associated with drug resistance.")
      
            
      val resistances = markers.filterNot(_._2.hasSusceptibility).filter(_._2.hasResistance).map(_._2.filter(_.isPresent))
      println(resistances.size + " SNPS detected associated with drug resistance.")
      resistances.foreach(println)
      
      //val susceptibilities = markers.filterNot(_.hasResistance).filter(_.hasSusceptibility).map(_.filter(_.isPresent))
      //susceptibilities.foreach(println)
      //println(susceptibilities.size)
      
      val nonDetected = markers.filterNot(_._2.hasResistance).filterNot(_._2.hasSusceptibility)
      //println(nonDetected.size + " non-detected genome positions with possible resistance.")
      //nonDetected.foreach(println)
      
      def ranges(ls: List[Int], range: Int): List[List[Int]] = ls match {
        case head :: tail => {
          val (b, a) = ls.span(p => p - head < 11)
          b:: ranges(a, 10)
        }
        case Nil => Nil
      }
      
      val within10bp = ranges(nonDetected.map(_._1).toList.sorted, 10).map(_.flatMap(nonDetected(_))).map(_.filter(_.markerType == "susceptibility"))
      println(within10bp.flatten.size + " non-detected SNP positions." )
      within10bp.foreach{ pList =>
        val drugs = pList.flatMap(_.drugs.split("_").toList).distinct
        val coordinates = pList.map(_.coordinate).mkString(", ")
        println("Possible resistance to " + drugs + " caused by one or more SNPs at genome coordinates " + coordinates + ".")
      }
      
      /**val within10bp = nonDetectedPos.groupBy{ c =>
          val idx = nonDetectedPos.indexOf(c)
          if (idx == 0) (nonDetectedPos(idx + 1) - nonDetectedPos(idx) < 11)
          else if (idx == nonDetectedPos.size - 1) (nonDetectedPos(idx) - nonDetectedPos(idx - 1) < 11)
          else (nonDetectedPos(idx + 1) - nonDetectedPos(idx) < 11) || (nonDetectedPos(idx) - nonDetectedPos(idx - 1) < 11)
      }
      println(within10bp.flatMap(_._2).size + " undetected SNP positions with possible resistance.")
      within10bp.foreach(println)
      println("Non-detected not within 10 bp.")
      within10bp(false).map(nonDetected.groupBy(_(0).coordinate)(_)).foreach(println)
      println("Non-detected within 10 bp.")
      val within10 = within10bp(true).map(nonDetected.groupBy(_(0).coordinate)(_))
      within10.foreach(println)*/
      
      val bothDetected = markers.filter(_._2.hasResistance).filter(_._2.hasSusceptibility)
      println(bothDetected.size + " genome positions with both resistance and susceptibility detected.")
      bothDetected.foreach(println)
      
      
      
      
    }
    
  }
}