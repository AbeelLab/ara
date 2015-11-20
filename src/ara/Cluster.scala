package ara

object Cluster {
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
  }
  
  
  
  
  implicit def seqtoBool(s: String) = new Cluster(s)
}