package ara

object Cluster extends MTBCclusters {
  class Cluster(s: String) {

    def hasReference: Boolean = s match {
      case "L1-L2-L3-L4" => true
      case "L2-L3-L4" => true
      case "L4" => true
      case "L4.2" => true
      case "L4.2.2" => true
      case "L4.2.2.2" => true
      case "L4.2.2.2.2" => true
      case "L4.2.2.2.2.2" => true
      case "L4.2.2.2.2.2.2" => true
      case "L4.2.2.2.2.2.2.2" => true
      case "L4.2.2.2.2.2.2.2.2" => true
      case "L4.2.2.2.2.2.2.2.2.2" => true
      case "L4.2.2.2.2.2.2.2.2.2.2" => true
      case _ => false
    }

    def getAncestor: String = s match {
      case "L1-L2-L3-L4" => "MTBC"
      case "L5-L6-LB" => "MTBC"
      case "L1" => "L1-L2-L3-L4"
      case "L2-L3-L4" => "L1-L2-L3-L4"
      case "L2-L3" => "L2-L3-L4"
      case "L4" => "L2-L3-L4"
      case "L2" => "L2-L3"
      case "L3" => "L2-L3"
      case "L5" => "L5-L6-LB"
      case "L6-LB" => "L5-L6-LB"
      case "L6" => "L6-LB"
      case "LB" => "L6-LB"
      case x if (x.contains(".")) => x.dropRight(2)
      case _ => "not known"
    }

    def hasAncestor: Boolean = s match {
      case "MTBC" => false
      case _ => true
    }
    
    def getSibling: String = s match {
      case x if (x.contains(".")) => if (x.endsWith("1")) x.dropRight(1) + "2" else x.dropRight(1) + "1"
      case "L1-L2-L3-L4" => "L5-L6-LB"
      case "L5-L6-LB" => "L1-L2-L3-L4"
      case "L1" => "L2-L3-L4"
      case "L2-L3-L4" => "L1"
      case "L2-L3" => "L4"
      case "L4" => "L2-L3"
      case "L2" => "L3"
      case "L3" => "L2"
      case "L5" => "L6-LB"
      case "L6-LB" => "L5"
      case "L6" => "LB"
      case "LB" => "L6"
      case _ => null
    }    
     
    def children: Array[String] = s match {
      case "MTBC" => Array("L1-L2-L3-L4", "L5-L6-LB")
      case "L1-L2-L3-L4" => Array("L1", "L2-L3-L4")
      case "L5-L6-LB" => Array("L5", "L6-LB")
      case "L2-L3-L4" => Array("L2-L3", "L4")
      case "L2-L3" => Array("L2", "L3")
      case "L6-LB" => Array("L6", "LB")
      case x if (mtbcClusters.contains(x)) => Array(x + ".1", x + ".2")
      case _ => Array(null, null)
    }
    
    def isCluster: Boolean = {
      mtbcClusters.contains(s)
    }
    
    def hasZeroMarkers: Boolean = {
      cNumbers(s)("markers") == 0
    }   

    def getSize: Int = {
      cNumbers(s)("size")
    }
    
  }

  implicit def seqtoBool(s: String) = new Cluster(s)
}